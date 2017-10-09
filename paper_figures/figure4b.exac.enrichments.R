#!/usr/bin/env Rscript

rm(list = ls())

dir = Sys.getenv('RAREVARDIR')
exacdir = Sys.getenv('EXAC_DIR')

# Load required packages
require(ggplot2)
require(gtable)
require(gridExtra)
require(RColorBrewer)
require(reshape2)
require(broom)
require(plyr)
require(dplyr)

#-------------------- FUNCTIONS

#-------------------- MAIN

#---- Global constants
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

# Define gene types to keep
types.to.keep = c('lincRNA', 'protein_coding')

# Outlier threshold 
medz.threshold = 2

#---- Analysis
# Read in table of ENSG and HGNC ids from GTEx V6P GTF
ensg.hgnc = read.table(paste0(dir, '/data/gtex.v6p.hgnc.ensg.mapper.txt'), 
  sep = '\t', header = F, stringsAsFactors = F)
names(ensg.hgnc) = c('gene', 'hgnc')

# Read in list of genes with type info
# Restrict to protein coding and lincRNA genes
ensg.types = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), 
  sep = '\t', header = F, stringsAsFactors = F)
colnames(ensg.types) = c('gene', 'type')
gene.info = ensg.types %>% 
  filter(type %in% types.to.keep) %>%
  merge(., ensg.hgnc, all.x = T)

# Read in ExAC file
# Calculate percentiles for each constraint score
# Restrict to genes in the intersection with our Gencode annotation
exac = read.table('/mnt/lab_data/montgomery/shared/ExAC/release0.3/functional_gene_constraint/forweb_cleaned_exac_r03_march16_z_data_pLI.txt', 
  sep = '\t', stringsAsFactors = F, header = T) %>%
  filter(gene %in% gene.info$hgnc)
names(exac)[2] = 'hgnc'

syn.quants = quantile(exac$syn_z, probs = seq(0, 1, .05))
mis.quants = quantile(exac$mis_z, probs = seq(0, 1, .05))
pli.quants = quantile(exac$pLI, probs = seq(0, 1, .1))

# Read in outliers data
# Add HGNC symbols 
medz = read.table(paste0(dir, '/data/outliers_medz_nothreshold_picked.txt'),
  sep = '\t', header = T, stringsAsFactors = F)
names(medz) = c('gene', 'indiv', 'ntissue', 'medz')
medz = medz %>% merge(., gene.info[, c('gene', 'hgnc')], all.x = T) %>%
  mutate(type = factor(gene.types[1], levels = gene.types))

# Read in Metatissue summary data
# Restrict to autosomal lincRNA and PC genes
# Add HGNC symbols to the Metatissue data
meta = read.table(paste0(dir, '/data/GTExReleaseV6PMetatissue.summary.txt'), 
  sep = '\t', header = T, stringsAsFactors = F)
names(meta) = tolower(names(meta))
meta = meta %>% filter(gene %in% gene.info$gene) %>%
  merge(., gene.info[, c('gene', 'hgnc')], all.x = T) %>%
  mutate(type = factor(gene.types[2], levels = gene.types))

# Read in GWAS genes data
gwas = read.delim(paste0(dir, '/features/annotations/GWAS/gwas_catalog_v1.0-downloaded_2015-11-30.tsv'), 
  header = T, sep = '\t', stringsAsFactors = F)

# Get list of reported genes
reported.genes = unique(sort(unlist(strsplit(gwas$REPORTED.GENE.S., ' - |, '))))

# Assign GWAS genes to the ExAC dataset
reported = exac %>% select(1,2,17,18,20) %>%
  mutate(annot = ifelse(hgnc %in% reported.genes, states[5], states[1]), type = factor(gene.types[4], levels = gene.types)) %>%
  select(hgnc, syn_z, mis_z, pLI, annot, type)

# Read in OMIM genes
omim.genes = read.table(paste0(dir, '/features/annotations/OMIM/omim.genes.txt'), 
  sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Assign OMIM genes to the ExAC dataset
omim = exac %>% select(1,2,17,18,20) %>%
  mutate(annot = ifelse(hgnc %in% omim.genes, states[6], states[1]), type = factor(gene.types[5], levels = gene.types)) %>%
  select(hgnc, syn_z, mis_z, pLI, annot, type)

# Add Synonymous and Missense Z-scores and pLI scores to the outlier and Metasoft datasets
# Remove genes without ExAC annotation
medz = merge(medz, exac[, c('hgnc', score.names)], all.x = T) %>%
  filter(!is.na(pLI))
meta = merge(meta, exac[, c('hgnc', score.names)], all.x = T) %>%
  filter(!is.na(pLI))

# Define outlier genes using magnitude threshold of 2
medz = medz %>% mutate(annot = ifelse(abs(medz) >= medz.threshold, states[2], states[1])) %>%
  select(hgnc, syn_z, mis_z, pLI, annot, type)

# Define top genes by P-value for RE2 model as shared eQTLs
# Take same number as the number of outlier genes
Ntop = table(medz$annot)[2]
meta = meta[order(meta$pvalue), ]
meta$annot = states[1]
meta$annot[1:Ntop] = states[3]
meta = meta %>% select(hgnc, syn_z, mis_z, pLI, annot, type)

# Analyze enrichment for loss-of-function, missense, and synonymous intolerant variation
# Compute odds ratios using Fisher's exact for GWAS, OMIM, eGenes, and outlier genes
syn.threshold = syn.quants['90%']
mis.threshold = mis.quants['90%']
pli.threshold = 0.9

binary.score.names = c('syn', 'mis', 'pli')

odds = rbind(medz, meta, reported, omim) %>%
  mutate(syn = syn_z > syn.threshold, mis = mis_z > mis.threshold, pli = pLI > pli.threshold) %>%
  melt(id.vars = c('hgnc', 'type', 'annot'), measure.vars = c('syn', 'mis', 'pli'), variable.name = 'score', value.name = 'annot.exac') %>%
  group_by(type, score) %>%
  do(tidy(fisher.test(table(.$annot, .$annot.exac))))

# Clean up the odds df for easy plotting
score.names.pretty = c('Synonymous intolerance', 'Missense intolerance', 'LoF intolerance')

odds = odds %>% ungroup() %>% mutate(score = mapvalues(score, from = binary.score.names, to = score.names.pretty), 
  type = factor(as.character(type), levels = gene.types[c(4, 5, 1, 2)]))

# Make the plot for pLI 
pli.odds.plot = ggplot(data = filter(odds, score == 'LoF intolerance'), aes(x = type, y = estimate)) + 
  geom_hline(yintercept = 1, colour = 'darkgrey') +
	geom_pointrange(aes(x = type, ymin = conf.low, ymax = conf.high, colour = type), size = .3) + theme_bw() +
	scale_colour_manual(values = type.colors) + xlab('') + ylab('Odds ratio') + guides(colour = F) +
	coord_flip()
pli.odds.plot

# Make the plot for other ExAC scores (supplemental)
fsize = 9

syn.mis.odds.plot = ggplot(data = filter(odds, score %in% score.names.pretty[1:2]), aes(x = type, y = estimate)) + 
  geom_hline(yintercept = 1, colour = 'darkgrey') +
  geom_pointrange(aes(x = type, ymin = conf.low, ymax = conf.high, colour = type), size = .3) + theme_bw() +
  scale_colour_manual(values = type.colors) + xlab('') + ylab('Odds ratio') + guides(colour = F) +
  coord_flip() + theme(strip.background = element_blank(),
          strip.text = element_text(size = fsize + 1, face = "bold"),
          axis.text = element_text(size = fsize),
          axis.title = element_text(size = fsize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(.5, 1, 1.5)) +
  facet_grid(. ~ score)
syn.mis.odds.plot

# TESTING
# Test for a significant difference in ExAC scores for eGenes and Outlier genes
# Merge the outlier and eGene dataset
# Create a predictor showing whether the gene is an eGene, an outlier gene, both or neither
# Run an ANOVA using this predictor for each score
# Perform Tukey's Range Test to see if there is a significant difference for just outlier genes compared to just eQTL genes
names(medz)[5] = 'outlier'
names(meta)[5] = 'egene'
out.egene.merged = merge(medz[, -6], meta[, -6])
out.egene.merged$combo = paste(out.egene.merged$outlier, out.egene.merged$egene, sep = ' - ')

# Synonymous Z-score
syn.test = TukeyHSD(aov(syn_z ~ combo, data = out.egene.merged))

# Missense Z-score
mis.test = TukeyHSD(aov(mis_z ~ combo, data = out.egene.merged))

# pLI
pli.test = TukeyHSD(aov(pLI ~ combo, data = out.egene.merged))

# Save workspace image
save.image(paste0(dir, '/data/figure4b.exac.enrichments.RData'))
