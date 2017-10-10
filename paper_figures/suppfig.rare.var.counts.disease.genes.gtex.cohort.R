#!/usr/bin/env Rscript

rm(list = ls())

dir = Sys.getenv('RAREVARDIR')
exacdir = Sys.getenv('EXAC_DIR')
annodir = paste0(dir, '/features/annotations')

# Load required packages
require(ggplot2)
require(reshape2)
require(plyr)
require(dplyr)
require(broom)
require(gtable)
require(ggbeeswarm)
require(gridExtra)
require(RColorBrewer)
require(cowplot)
require(data.table)

#------------- FUNCTIONS

#------------- MAIN

# Global constants
var.types = c('SNV', 'Indel', 'SV')

maf.bins = c(0, 0.01, 0.05, 0.1, 0.25)
maf.bins.labels = c('0-1%', '1-5%', '5-10%', '10-25%')
maf.bins.alphas = seq(1, .4, -.2)
names(maf.bins.alphas) = maf.bins.labels

types.to.keep = c('lincRNA', 'protein_coding')

medz.threshold = 2

out.states = c('Control', 'Outlier')

egene.states = c('Control', 'eGene')

cardio.states = c('Control', 'Cardio')

acmg.states = c('Control', 'ACMG')

gwas.states = c('Control', 'GWAS')

cancer.states = c('Control', 'Cancer')

orpha.states = c('Control', 'Orphanet')

omim.states = c('Control', 'OMIM')

clin.states = c('Control', 'ClinVar')

lof.states = c('Control', 'LoF intolerant')

ddg2p.states = c('Control', 'DDG2P')

gene.list.colours = data.frame(List = c('Control', 'Outlier', 'eGene', 'Cardio', 'ACMG', 'GWAS', 'Cancer', 'Orphanet', 'OMIM', 'ClinVar', 'LoF intolerant', 'DDG2P'), 
	Colour = c('darkgrey', 'dodgerblue3', '#A6CEE3', '#FF7F00', '#6A3D9A', '#33A02C', '#B15928', '#E31A1C', '#FB9A99', '#FDBF6F', '#FF7F00', 'springgreen4'), stringsAsFactors = F)
list.colours = gene.list.colours$Colour
names(list.colours) = gene.list.colours$List

fsize = 11
axisFontSizes = theme(axis.text.x = element_text(size = fsize),
                      axis.text.y = element_text(size = fsize),
                      axis.title = element_text(size = fsize),
                      strip.text = element_text(size = fsize),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())

# Read in the gene position, type, and HGNC ID info
gene.data = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs.bed'), 
	sep = '\t', header = F, stringsAsFactors = F)
names(gene.data) = c('chrom', 'start', 'stop', 'gene')

ensg.types = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), 
	sep = '\t', header = F, stringsAsFactors = F)
colnames(ensg.types) = c('gene', 'gene.type')
ensg.types = ensg.types %>% filter(gene.type %in% types.to.keep)

ensg.hgnc = read.table(paste0(dir, '/data/gtex.v6p.hgnc.ensg.mapper.txt'), sep = '\t', header = F, stringsAsFactors = F)
colnames(ensg.hgnc) = c('gene', 'hgnc')

gene.data = merge(ensg.types, gene.data, all.x = T) %>% merge(., ensg.hgnc, all.x = T) %>%
	select(chrom, start, stop, gene, hgnc, gene.type)

#-- Process the variant data

# Read in the variant data for MAFs 0-25%
# Will have gene and ER promoter annotation
vars = fread(paste0('zcat ', dir, '/features/variantBeds/gtex.vars.maf0-25.er.prom.bed.gz')) %>%
	setNames(c('chrom', 'start', 'stop', 'maf', 'gene', 'var.type', 'promoter')) %>%
	filter(gene %in% gene.data$gene) %>%
	mutate(var.type = mapvalues(var.type, from = c('SNPs', 'indels', 'HallLabSV'), to = var.types))

# Generate count info by variant type for each gene
# If the gene has no variant in a given MAF bin and variant type, fill in with 0
exp.design = expand.grid(gene.data$gene, maf.bins.labels)
names(exp.design) = c('gene', 'maf.bin')
exp.design = merge(exp.design, gene.data[, c('gene', 'hgnc', 'gene.type')], all.x = T)

var.counts = vars %>% 
	mutate(maf.bin = cut(maf, breaks = maf.bins, labels = maf.bins.labels)) %>%
	group_by(gene, maf.bin) %>%
	summarise(count = n()) %>%
	ungroup() %>%
	merge(., exp.design, all = T) %>%
	mutate(count = ifelse(is.na(count), 0, count), pretty = gsub('\\.[0-9]+', '', gene))

# Add gene set annotations
# Outlier genes
medz = read.table(paste0(dir, '/data/outliers_medz_nothreshold_picked.txt'), sep = '\t', header = T, stringsAsFactors = F)
names(medz) = c('gene', 'indiv', 'ntissue', 'medz')
outlier.genes = medz %>% mutate(Outlier = factor(out.states[1 + (abs(medz) >= medz.threshold)], levels = out.states)) %>%
	filter(Outlier == out.states[2]) %>% select(gene)
var.counts = var.counts %>% mutate(Outlier = factor(out.states[1 + (gene %in% outlier.genes$gene)], levels = out.states))

# Mutli-tissue eGenes
meta = read.table(paste0(dir, '/data/GTExReleaseV6PMetatissue.summary.txt'), sep = '\t', header = T, stringsAsFactors = F)
names(meta)[1] = 'gene'
meta = meta %>% filter(gene %in% var.counts$gene)
meta = meta[order(meta$Pvalue), ]
meta$egene = 0
meta$egene[1:nrow(outlier.genes)] = 1
meta.genes = meta %>% filter(egene == 1) %>% select(gene)
var.counts = var.counts %>% mutate(eGene = factor(egene.states[1 + (gene %in% meta.genes$gene)], levels = egene.states))

# Cardio
cardio.gold.genes = read.table(paste0(annodir, '/GeneListsOther/cardio.genes.gold.standard.csv'), sep = ',', header = T, stringsAsFactors = F)[, 1]
var.counts = var.counts %>% mutate(Cardio = factor(cardio.states[1 + (hgnc %in% cardio.gold.genes)], levels = cardio.states))

# ACMG
acmg = read.table(paste0(annodir, '/ACMG/completed_acmg.csv'), sep = ',', header = T, stringsAsFactors = F)
var.counts = var.counts %>% mutate(ACMG = factor(acmg.states[1 + (hgnc %in% acmg$Gene)], levels = acmg.states))

# GWAS
gwas = read.delim(paste0(annodir, '/GWAS/gwas_catalog_v1.0-downloaded_2015-11-30.tsv'), header = T, sep = '\t', stringsAsFactors = F)
reported.genes = unique(sort(unlist(strsplit(gwas$REPORTED.GENE.S., ' - |, '))))
var.counts = var.counts %>% mutate(GWAS = factor(gwas.states[1 + (hgnc %in% reported.genes)], levels = gwas.states))

# Cancer
cancer.gold.genes = read.table(paste0(annodir, '/GeneListsOther/cancer.genes.gold.standard.csv'), sep = ',', header = T, stringsAsFactors = F)[, 1]
var.counts = var.counts %>% mutate(Cancer = factor(cancer.states[1 + (hgnc %in% cancer.gold.genes)], levels = cancer.states))

# Orphanet
orpha.genes = read.table(paste0(annodir, '/OrphaNet/orphanet.genes.txt'), sep = '\t', header = F, stringsAsFactors = F)[, 1]
var.counts = var.counts %>% mutate(Orphanet = factor(orpha.states[1 + (pretty %in% orpha.genes)], levels = orpha.states))

# OMIM
omim.genes = read.table(paste0(annodir, '/OMIM/omim.genes.txt'), sep = '\t', header = F, stringsAsFactors = F)[, 1]
var.counts = var.counts %>% mutate(OMIM = factor(omim.states[1 + (hgnc %in% omim.genes)], levels = omim.states))

# ClinVar
clin = read.csv(paste0(annodir, '/GeneListsOther/gene_condition_source_id'), sep = '\t', stringsAsFactors = F)
names(clin)[2] = 'Gene'
var.counts = var.counts %>% mutate(ClinVar = factor(clin.states[1 + (hgnc %in% clin$Gene)], levels = clin.states))

# LoF intolerant genes from ExAC
# Read in and define LoF intolerant genes (pLI > 0.9)
exac.threshold = 0.9
lof.intolerant.genes = read.table(paste0(exacdir, '/forweb_cleaned_exac_r03_march16_z_data_pLI.txt'), 
	sep = '\t', stringsAsFactors = F, header = T) %>% 
	filter(pLI > exac.threshold) %>%
	select(gene)
var.counts = var.counts %>% mutate(`LoF intolerant` = factor(lof.states[1 + (hgnc %in% lof.intolerant.genes$gene)], levels = lof.states))

# DDG2P
ddg2p = read.table(paste0(annodir, '/DDG2P/DDG2P_2_8_2017.csv.gz'), sep = ',', header = T, stringsAsFactors = F)
ddg2p.genes = unique(ddg2p$gene.symbol)
var.counts = var.counts %>% mutate(DDG2P = factor(ddg2p.states[1 + (hgnc %in% ddg2p.genes)], levels = ddg2p.states))

#---- Test for differences using a t-test
# Compute difference in rare variant counts between disease and control genes for each variant type and each gene list
tstats = melt(var.counts, measure.vars = names(var.counts)[names(var.counts) %in% names(list.colours)], 
	variable.name = 'list.type', value.name = 'list') %>%
	mutate(list = factor(list, levels = rev(names(list.colours)))) %>%
	group_by(list.type, maf.bin) %>% 
	do(tidy(t.test(count ~ list, data = .))) %>%
	ungroup() %>% 
	mutate(list.type = factor(list.type, levels = names(list.colours)),
		maf.bin = factor(maf.bin, levels = maf.bins.labels))

# Make plot of the rare variant count differences
tstats.plot = ggplot(data = filter(tstats, list.type != 'LoF intolerant'), aes(x = list.type, y = estimate)) + theme_bw() + 
	geom_hline(yintercept = 0, colour = 'darkgrey') + geom_pointrange(aes(x = list.type, ymin = conf.low, ymax = conf.high, 
		colour = list.type, alpha = maf.bin), size = .3, position = position_dodge(width = .7)) + 
	theme(legend.position = c(.8, .8), strip.background = element_blank()) + xlab('') + 
	ylab('Difference in mean number of variants\nbetween annotated and control genes') + scale_colour_manual(values = list.colours) + 
	guides(colour = F) + 
	scale_alpha_manual(values = maf.bins.alphas, name = 'MAF') +
	axisFontSizes
tstats.plot

# Save workspace image
save.image(paste0(dir, '/data/suppfig.rare.var.counts.disease.genes.gtex.cohort.RData'))




