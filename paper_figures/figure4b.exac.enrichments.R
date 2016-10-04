#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')
exacdir = Sys.getenv('EXAC_DIR')

# Load required packages
require(ggplot2)
require(RColorBrewer)

#-------------------- FUNCTIONS

# Function to find Synonymous, Missense and pLI scores for a given gene
addScores <- function(x, exac){
  index = which(exac$gene == x)
  if(length(index) == 0){
    out = rep(NA, 3)
  }
  else{
    out = c(exac$syn_z[index], exac$mis_z[index], exac$pLI[index])
  }
  return(out)
}

# Quantile normalization from Mauro
getInverseNormal <- function(x){
        x <- as.numeric(x)
        xx2<-qnorm((rank(x,na.last = "keep") - 0.5) / sum(!is.na(x)))
        return(xx2)
}

#-------------------- MAIN

# Define gene types and colors
states = c('Background', 
  'Genes with an outlier', 
  'Genes with a shared eQTL', 
  'Mapped GWAS gene', 
  'Reported GWAS gene', 
  'OMIM gene', 
  'Orphanet gene')

colors = brewer.pal(12, name = 'Paired')
colors = c('darkgrey', 'dodgerblue3', colors[c(1, 3:6)])

names(colors) = states

gene.types = c('Outlier', 'eGene', 'GWAS mapped', 'GWAS', 'OMIM', 'Orphanet')
type.colors = colors[-1]
names(type.colors) = gene.types

# Define ExAC variable names and colors
score.names = c('syn_z', 'mis_z', 'pLI')
bin.vars = c('syn_z_quart', 'mis_z_quart', 'pLI_quart')
exac.colors = c('darkgrey', 'orangered3', 'orangered4')
names(exac.colors) = bin.vars

# Read in table of ENSG and HGNC ids from GTEx V6P GTF
ensg.hgnc = read.table(paste0(dir, '/data/gtex.v6p.hgnc.ensg.mapper.txt'), sep = '\t', header = F, stringsAsFactors = F, row.names = 1)

# Read in list of genes with type info
# Restrict to protein coding and lincRNA genes
ensg.types = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), sep = '\t', header = F, stringsAsFactors = F)
colnames(ensg.types) = c('GENE', 'TYPE')
ensg.types$HGNC = ensg.hgnc[ensg.types[, 1], 1]
rownames(ensg.types) = ensg.types$GENE
types.to.keep = c('lincRNA', 'protein_coding')
ensg.types = ensg.types[ensg.types$TYPE %in% types.to.keep, ]
ensg.types$PRETTY = gsub('\\.*', '', ensg.types$GENE)

# Read in ExAC file
# Calculate percentiles for each constraint score
# Restrict to genes in the intersection with our Gencode annotation
exac = read.table(paste0(exacdir, '/forweb_cleaned_exac_r03_march16_z_data_pLI.txt'), sep = '\t', stringsAsFactors = F, header = T)
exac = exac[exac$gene %in% ensg.types$HGNC, ]
syn.quants = quantile(exac$syn_z, probs = seq(0, 1, .05))
mis.quants = quantile(exac$mis_z, probs = seq(0, 1, .05))
pli.quants = quantile(exac$pLI, probs = seq(0, 1, .1))

# Read in outliers data
medz = read.table(paste0(dir, '/data/outliers_medz_nothreshold_picked.txt'), sep = '\t', header = T, stringsAsFactors = F)

# Add HGNC symbols to medz dataset
medz$HGNC = unlist(lapply(medz$GENE, function(x) ensg.types[x, 'HGNC']))

# Read in Metasoft summary data
# Restrict to autosomal lincRNA and PC genes
meta = read.table(paste0(dir, '/data/GTExReleaseV6PMetasoft.summary.txt'), sep = '\t', header = T, stringsAsFactors = F)
meta = meta[meta$Gene %in% ensg.types$GENE, ]

# Add HGNC symbols to the metasoft dataset
meta$HGNC = unlist(lapply(meta$Gene, function(x) ensg.types[x, 'HGNC']))

# Read in GWAS genes data
gwas = read.delim(paste0(dir, '/features/annotations/GWAS/gwas_catalog_v1.0-downloaded_2015-11-30.tsv'), header = T, sep = '\t', stringsAsFactors = F)

# Get list of reported genes
reported.genes = unique(sort(unlist(strsplit(gwas$REPORTED.GENE.S., ' - |, '))))

# Assign GWAS genes to the ExAC dataset
reported = exac[, c(1,2,17,18,20)]
reported$GWAS = ifelse(exac$gene %in% reported.genes, states[5], states[1])

# Read in OMIM genes
omim.genes = read.table(paste0(dir, '/features/annotations/OMIM/omim.genes.txt'), sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Assign OMIM genes to the ExAC dataset
omim = exac[, c(1,2,17,18,20)]
omim$OMIM = ifelse(exac$gene %in% omim.genes, states[6], states[1])

# Add Synonymous and Missense Z-scores and pLI scores to the outlier and Metasoft datasets
medz[, score.names] = matrix(unlist(lapply(medz$HGNC, addScores, exac = exac)), ncol = 3, byrow = T)
meta[, score.names] = matrix(unlist(lapply(meta$HGNC, addScores, exac = exac)), ncol = 3, byrow = T)

# Remove genes without ExAC annotation
medz = medz[!is.na(medz$pLI), ]
meta = meta[!is.na(meta$pLI), ]

# Define outlier genes using magnitude threshold of 2
medz.threshold = 2
medz$OUTLIER = states[1]
medz$OUTLIER[abs(medz$Z) >= medz.threshold] = states[2]

# Define top genes by P-value for RE2 model as shared eQTLs
# Take same number as the number of outlier genes
Ntop = table(medz$OUTLIER)[2]
meta = meta[order(meta$Pvalue), ]
meta$EQTL = states[1]
meta$EQTL[1:Ntop] = states[3]

# Analyze enrichment for PTV, missense, and synonymous intolerant variation
# Restrict analysis to PTV-intolerant genes 
syn.threshold = syn.quants['90%']
mis.threshold = mis.quants['90%']
pli.threshold = 0.9

medz = medz[, c('HGNC', score.names, 'OUTLIER')]
medz$SYN = medz$syn_z > syn.threshold
medz$MIS = medz$mis_z > mis.threshold
medz$PTV = medz$pLI > pli.threshold

meta = meta[, c('HGNC', score.names, 'EQTL')]
meta$SYN = meta$syn_z > syn.threshold
meta$MIS = meta$mis_z > mis.threshold
meta$PTV = meta$pLI > pli.threshold

reported = reported[, c('gene', score.names, 'GWAS')]
names(reported)[1] = 'HGNC'
reported$SYN = reported$syn_z > syn.threshold
reported$MIS = reported$mis_z > mis.threshold
reported$PTV = reported$pLI > pli.threshold

omim = omim[, c('gene', score.names, 'OMIM')]
names(omim)[1] = 'HGNC'
omim$SYN = omim$syn_z > syn.threshold
omim$MIS = omim$mis_z > mis.threshold
omim$PTV = omim$pLI > pli.threshold

# Compute odds ratios using Fisher's exact for GWAS, OMIM, eGenes, and outlier genes
gwas.odds.syn = fisher.test(table(reported$GWAS, reported$SYN))
gwas.odds.mis = fisher.test(table(reported$GWAS, reported$MIS))
gwas.odds.ptv = fisher.test(table(reported$GWAS, reported$PTV))

omim.odds.syn = fisher.test(table(omim$OMIM, omim$SYN))
omim.odds.mis = fisher.test(table(omim$OMIM, omim$MIS))
omim.odds.ptv = fisher.test(table(omim$OMIM, omim$PTV))

egenes.odds.syn = fisher.test(table(meta$EQTL, meta$SYN))
egenes.odds.mis = fisher.test(table(meta$EQTL, meta$MIS))
egenes.odds.ptv = fisher.test(table(meta$EQTL, meta$PTV))

outlier.odds.syn = fisher.test(table(medz$OUTLIER, medz$SYN))
outlier.odds.mis = fisher.test(table(medz$OUTLIER, medz$MIS))
outlier.odds.ptv = fisher.test(table(medz$OUTLIER, medz$PTV))

# Make data frames for easy plotting
score.names.pretty = c('Synonymous intolerance', 'Missense intolerance', 'PTV intolerance')

ptv.odds = data.frame(SET = factor(gene.types[c(4, 5, 1, 2)], levels = gene.types[c(4, 5, 1, 2)]), 
			POINT = c(gwas.odds.ptv$estimate, omim.odds.ptv$estimate, outlier.odds.ptv$estimate, egenes.odds.ptv$estimate), 
      PVAL = c(gwas.odds.ptv$p.value, omim.odds.ptv$p.value, outlier.odds.ptv$p.value, egenes.odds.ptv$p.value))
ptv.odds[, c('LOW', 'HIGH')] = matrix(c(gwas.odds.ptv$conf.int, omim.odds.ptv$conf.int, outlier.odds.ptv$conf.int, egenes.odds.ptv$conf.int), ncol = 2, byrow = T)
ptv.order = order(ptv.odds$POINT)
ptv.odds = ptv.odds[ptv.order, ]
ptv.odds$SET = factor(as.character(ptv.odds$SET), levels = rev(as.character(ptv.odds$SET)))
ptv.odds$SCORE = factor(score.names.pretty[3], levels = score.names.pretty)

syn.odds = data.frame(SET = factor(gene.types[c(4, 5, 1, 2)], levels = gene.types[c(4, 5, 1, 2)]), 
      POINT = c(gwas.odds.syn$estimate, omim.odds.syn$estimate, outlier.odds.syn$estimate, egenes.odds.syn$estimate), 
      PVAL = c(gwas.odds.syn$p.value, omim.odds.syn$p.value, outlier.odds.syn$p.value, egenes.odds.syn$p.value))
syn.odds[, c('LOW', 'HIGH')] = matrix(c(gwas.odds.syn$conf.int, omim.odds.syn$conf.int, outlier.odds.syn$conf.int, egenes.odds.syn$conf.int), ncol = 2, byrow = T)
syn.odds = syn.odds[ptv.order, ]
syn.odds$SET = factor(as.character(syn.odds$SET), levels = rev(as.character(syn.odds$SET)))
syn.odds$SCORE = factor(score.names.pretty[1], levels = score.names.pretty)

mis.odds = data.frame(SET = factor(gene.types[c(4, 5, 1, 2)], levels = gene.types[c(4, 5, 1, 2)]), 
      POINT = c(gwas.odds.mis$estimate, omim.odds.mis$estimate, outlier.odds.mis$estimate, egenes.odds.mis$estimate), 
      PVAL = c(gwas.odds.mis$p.value, omim.odds.mis$p.value, outlier.odds.mis$p.value, egenes.odds.mis$p.value))
mis.odds[, c('LOW', 'HIGH')] = matrix(c(gwas.odds.mis$conf.int, omim.odds.mis$conf.int, outlier.odds.mis$conf.int, egenes.odds.mis$conf.int), ncol = 2, byrow = T)
mis.odds = mis.odds[ptv.order, ]
mis.odds$SET = factor(as.character(mis.odds$SET), levels = rev(as.character(mis.odds$SET)))
mis.odds$SCORE = factor(score.names.pretty[2], levels = score.names.pretty)

full.odds = rbind(ptv.odds, syn.odds, mis.odds)

# Make the plot for PTV intolerance
ptv.odds.plot = ggplot(data = ptv.odds, aes(x = SET, y = POINT)) + geom_hline(yintercept = 1, colour = 'darkgrey') +
	geom_pointrange(aes(x = SET, ymin = LOW, ymax = HIGH, colour = SET), size = .3) + theme_bw() +
	scale_colour_manual(values = type.colors) + xlab('') + ylab('Odds ratio') + guides(colour = F) +
	coord_flip() #+ ggtitle('PTV-intolerance')

# Make the plot for other ExAC scores (supplemental)
pdf(paste0(dir, '/paper_figures/suppfig.exac.enrichments.other.scores.pdf'), height = 2.5, width = 4.5)

fsize = 9

full.odds.plot = ggplot(data = full.odds[full.odds$SCORE %in% score.names.pretty[1:2], ], aes(x = SET, y = POINT)) + geom_hline(yintercept = 1, colour = 'darkgrey') +
  geom_pointrange(aes(x = SET, ymin = LOW, ymax = HIGH, colour = SET), size = .3) + theme_bw() +
  scale_colour_manual(values = type.colors) + xlab('') + ylab('Odds ratio') + guides(colour = F) +
  coord_flip() + theme(strip.background = element_blank(),
          strip.text = element_text(size = fsize + 1, face = "bold"),
          axis.text = element_text(size = fsize),
          axis.title = element_text(size = fsize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(.5, 1, 1.5)) +
  facet_grid(. ~ SCORE)
full.odds.plot

dev.off()

# TESTING
# Test for a significant difference in ExAC scores for eGenes and Outlier genes
# Merge the outlier and eGene dataset
# Create a predictor showing whether the gene is an eGene, an outlier gene, both or neither
# Run an ANOVA using this predictor for each score
# Perform Tukey's Range Test to see if there is a significant difference for just outlier genes compared to just eQTL genes
out.egene.merged = merge(medz, meta)
out.egene.merged$COMBO = paste(out.egene.merged$OUTLIER, out.egene.merged$EQTL, sep = ' - ')

# Synonymous Z-score
syn.test = TukeyHSD(aov(syn_z ~ COMBO, data = out.egene.merged))

# Missense Z-score
mis.test = TukeyHSD(aov(mis_z ~ COMBO, data = out.egene.merged))

# pLI
pli.test = TukeyHSD(aov(pLI ~ COMBO, data = out.egene.merged))

# Save workspace image
save.image(paste0(dir, '/data/figure4b.exac.enrichments.RData'))
