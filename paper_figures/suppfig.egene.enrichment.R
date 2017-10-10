#!/usr/bin/env Rscript

rm(list = ls())

dir = Sys.getenv('RAREVARDIR')
gtexcisdir = Sys.getenv('GTEXCISDIR')

# Load required packages
require(ggplot2)
require(gtable)
require(gridExtra)
require(cowplot)
require(dplyr)
#-------------------- FUNCTIONS



#-------------------- MAIN

# Define outlier states and colors
# Define colors for the ExAC comparison
out.states = c('Background', 'Genes with an outlier')

gtex.colors = c('darkgrey', 'dodgerblue3')
names(gtex.colors) = out.states

# Define quantile names
decile.names = paste(seq(10, 100, 10))

# Read in list of genes with type info
# Restrict to protein coding and lincRNA genes
ensg.types = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), 
  sep = '\t', header = F, stringsAsFactors = F)
colnames(ensg.types) = c('GENE', 'TYPE')
types.to.keep = c('lincRNA', 'protein_coding')
ensg.types = filter(ensg.types, TYPE %in% types.to.keep) %>%
  mutate(PRETTY = gsub('\\.[0-9]+', '', GENE))

# Read in most extreme outliers for each gene
medz = read.table(paste0(dir, '/data/outliers_medz_nothreshold_picked.txt'), 
  sep = '\t', header = T, stringsAsFactors = F)

# Define outlier genes 
medz.threshold = 2
outlier.genes = filter(medz, abs(Z) >= medz.threshold) %>%
  select(GENE)

# Read in list of tissues
tissues = read.table(paste0(dir, 'preprocessing/gtex_2015-01-12_tissues_all_normalized_samples.txt'), 
  sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Read in median expression across tissues
medians = read.table(paste0(dir, '/preprocessing/gtex_2015-01-12_median_rpkm.txt'), 
  header = T, stringsAsFactors = F)
medians = filter(medians, Gene %in% ensg.types$GENE)
rownames(medians) = medians$Gene

# Make matrix of genes by tissues to hold eGene occurrences
# Will use an FDR threshold of 0.05
FDR = 0.05
egenes = matrix(0, ncol = length(tissues), nrow = nrow(medz))
rownames(egenes) = medz$GENE
colnames(egenes) = tissues

for(i in 1:length(tissues)){
  tissue_egenes = read.table(paste(gtexcisdir, tissues[i], '_Analysis.v6p.FOR_QC_ONLY.egenes.txt.gz', sep = ''), 
    sep = '\t', header = T, stringsAsFactors = F)
  tissue_egenes = filter(tissue_egenes, gene_id %in% ensg.types$GENE) %>%
    filter(qval <= FDR) %>%
    filter(gene_id %in% rownames(egenes))
  egenes[tissue_egenes$gene_id, tissues[i]] = 1
}

# Make dataframe of genes with number of tissues each gene is an eGene in along with outlier status and mean RPKM across tissues
egenes.df = data.frame(GENE = rownames(egenes), COUNT = rowSums(egenes), stringsAsFactors = F)
egenes.df = mutate(egenes.df, OUTLIER = factor(ifelse(GENE %in% outlier.genes[, 1], out.states[2], out.states[1]), levels = out.states)) %>%
  mutate(RPKM = medians[GENE, 'Median'])

# Add info about mean RPKM decile
gtex.deciles = quantile(medians$Median, prob = seq(0, 1, .1), na.rm = T)
egenes.df = mutate(egenes.df, DECILE = cut(egenes.df$RPKM, breaks = gtex.deciles, include.lowest = T, labels = decile.names))

# Read in Metatissue results
meta = read.table(paste0(dir, '/data/GTExReleaseV6PMetatissue.summary.txt'), 
  sep = '\t', header = T, stringsAsFactors = F)
meta = meta %>% select(Gene, eGenes) %>% filter(Gene %in% ensg.types$GENE)
colnames(meta) = c('GENE', 'COUNT')
rownames(meta) = meta$GENE

# Add outlier status and RPKM info to Metatissue data
meta = meta %>% mutate(OUTLIER = factor(ifelse(GENE %in% outlier.genes[, 1], out.states[2], out.states[1]), levels = out.states)) %>%
  mutate(RPKM = medians[GENE, 'Median'], DECILE = cut(RPKM, breaks = gtex.deciles, include.lowest = T, labels = decile.names))

# Remove genes without RPKM info
meta = filter(meta, !is.na(RPKM))
egenes.df = filter(egenes.df, !is.na(RPKM))

# Combine Metasoft and TBT results
egenes.df$TYPE = 'Tissue by tissue'
meta$TYPE = 'Meta-Tissue'

full = rbind(egenes.df, meta)
full$TYPE = factor(full$TYPE, levels = c('Meta-Tissue', 'Tissue by tissue'))

# Make density plot version of above
egenes.meta.dens = ggplot(data = full, aes(x = COUNT, fill = OUTLIER)) + 
    geom_density(alpha = .8) + theme_bw() + xlab('Number of tissues with an eQTL') + ylab('Density') + 
    scale_fill_manual(values = gtex.colors) + ggtitle('GTEx eGenes') + 
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.position = 'top',
          strip.background = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10), 
          legend.text = element_text(size = 10),
          legend.key = element_blank()) + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = .5)) + facet_grid(. ~ TYPE)
egenes.meta.dens

# Make supplementary plot showing effect of mean RPKM across tissues on eGene count distribution stratified by outlier status
# For TBT and Metasoft
egenes.meta.rpkm.box = ggplot(data = meta, aes(x = DECILE, y = COUNT, fill = OUTLIER)) + 
  geom_boxplot(outlier.colour = NA) + theme_bw() + xlab('Mean RPKM bin (%)') +
  ylab('Number of tissues') + scale_fill_manual(values = gtex.colors) + 
  ggtitle('GTEx eGenes') +
  guides(fill = guide_legend(title = NULL)) + theme(legend.position = 'top', 
    axis.text = element_text(size = 8), axis.title = element_text(size = 10), 
    legend.text = element_text(size = 10), legend.key = element_blank()) + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = .5))
egenes.meta.rpkm.box

# Test for significant difference between genes w/ and w/o outliers in terms of eGene count across tissues
# Will use Wilcoxon Rank Sum test
# Also calculate median values for background and outliers
egenes.test = wilcox.test(COUNT ~ OUTLIER, data = egenes.df)
egenes.test

egenes.med.back = summary(egenes.df$COUNT[egenes.df$OUTLIER == out.states[1]])
egenes.med.back

egenes.med.out = summary(egenes.df$COUNT[egenes.df$OUTLIER == out.states[2]])
egenes.med.out

meta.test = wilcox.test(COUNT ~ OUTLIER, data = meta)
meta.test

meta.med.back = summary(meta$COUNT[meta$OUTLIER == out.states[1]])
meta.med.back

meta.med.out = summary(meta$COUNT[meta$OUTLIER == out.states[2]])
meta.med.out

# Repeat same test as above but accounting for RPKM decile
meta.test.rpkm.pvals = c()
for(i in 1:length(decile.names)){
  meta.test.rpkm.pvals[i] = wilcox.test(COUNT ~ OUTLIER, data = meta[meta$DECILE == decile.names[i], ])$p.value
}

# make combined supplementary figure with the exac scores and disease genes
load(paste0(dir, '/data/figure4b.exac.enrichments.RData'))
load(paste0(dir, '/data/figure4c.gene.list.enrichments.RData'))
load(paste0(dir, '/data/suppfig.rare.var.counts.disease.genes.gtex.cohort.RData'))

fsize = 9
fontSizes = theme(legend.key.height = unit(0.8, "line"),
                  legend.key = element_blank(),
                  legend.text = element_text(size = fsize-1),
                  legend.background = element_blank(),
                  axis.text.x = element_text(size = fsize),
                  axis.text.y = element_text(size = fsize),
                  axis.title.x = element_text(size = fsize),
                  axis.title.y = element_text(size = fsize),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = fsize + 1),
                  strip.text = element_text(size = fsize, face = "bold"))

exac.supp.plot = syn.mis.odds.plot + ggtitle('ExAC constraint') + theme(plot.title = element_text(hjust = .5)) + fontSizes

full.ors.plot.alt = full.ors.plot.alt + fontSizes + theme(legend.position = c(0.75, 0.85),
                                                  plot.margin = margin(5.5,5.5,5.5,-5),
                                                  legend.key.width = unit(0.6, "line"), 
                                                  legend.title = element_blank())

egenes.meta.dens = egenes.meta.dens + fontSizes + theme(legend.position = c(0.25, 0.8)) + theme(legend.title = element_blank())

egenes.meta.rpkm.box = egenes.meta.rpkm.box + ggtitle('') + fontSizes + guides(fill = FALSE) +
    ylab('Number of tissues with an eQTL')

tstats.plot = tstats.plot + theme(legend.position = c(.8, .7), legend.title = element_text(size = fsize)) + fontSizes

suppfig.evol.constraint.egene.enrich.combined = ggdraw() +
    draw_plot(exac.supp.plot, 0, 2/3, 0.6, 1/3) +
    draw_plot(full.ors.plot.alt, 0.6, 2/3, 0.4, 1/3) +
    draw_plot(tstats.plot, 0.05, 1/3, .95, 1/3) +
    draw_plot(egenes.meta.dens, 0, 0, 0.6, 1/3) +
    draw_plot(egenes.meta.rpkm.box, 0.6, 0, 0.4, 1/3 - .05) +
    draw_plot_label(c('a','b','c','d','e'), c(0,0.6,0,0,0.6), c(1,1,2/3,1/3,1/3), size = 11)

pdf(paste0(dir, '/paper_figures/suppfig.exac.egenes.egenesrpkm.pdf'), height = 7.5, width = 7.2)

suppfig.evol.constraint.egene.enrich.combined

dev.off()

# Save workspace image
save.image(paste0(dir, '/data/suppfig.egene.enrichments.RData'))

