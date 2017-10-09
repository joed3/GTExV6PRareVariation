#!/usr/bin/env Rscript
rm(list = ls())

dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(ggplot2)
require(gtable)
require(gridExtra)
require(dplyr)
require(plyr)
require(reshape2)
require(bedr)

#----------------- FUNCTIONS


#----------------- MAIN

# Constants
out.states = c('Non-outlier', 'Outlier')
out.states.strat = c('Non-outlier', 'Under', 'Over')

out.colors = c('darkgrey', 'dodgerblue3')
names(out.colors) = out.states

strat.colors = c('darkgrey', 'dodgerblue1', 'dodgerblue3')
names(strat.colors) = out.states.strat

# Read in list of individuals with WGS data
wgs.inds = read.table(paste0(dir, '/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt'), sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Read in list of outliers
outliers = read.table(paste0(dir, '/shared_data/outliers_medz_picked.txt'), sep = '\t', header = T, stringsAsFactors = F) %>%
	mutate(Outlier = out.states[2])

# Subset outliers to those with WGS data
outliers = outliers %>% filter(INDS %in% wgs.inds)

# Read in full set of Z-scores and subset to the genes with outliers
# Restrict to the WGS individuals with data for each gene
medz = read.table(paste0(dir, '/shared_data/outliers_medz_zscores.txt'), sep = '\t', header = T, stringsAsFactors = F)
colnames(medz) = gsub('\\.', '-', colnames(medz))
medz = medz %>% filter(GENE %in% outliers$GENE) %>% select(which(colnames(medz) %in% c('GENE', wgs.inds))) %>%
	melt(data = ., id.vars = 'GENE', value.name = 'Z', variable.name = 'INDS') %>% 
	filter(!is.na(Z)) %>% merge(., outliers[, c(1:2, 5)], all.x = T) %>%
	mutate(Outlier = ifelse(is.na(Outlier), out.states[1], Outlier))
colnames(medz)[1:2] = c('Gene', 'Ind')

# Read in list of all rare variants (MAF < 1%) near each gene (<= 10kb from the TSS/TES) in each individual
var.data = read.table(paste0(dir, '/features/byGene/10kb/all_rare_variants_SNPs.txt'), sep = '\t', header = F, stringsAsFactors = F)
colnames(var.data) = c('Ind', 'Gene', 'chr', 'pos')

# Merge the outlier data with the rare variant data
var.data = merge(medz, var.data, by = c('Gene', 'Ind'))

# Add promoter annotations from Epigenomics Roadmap
# Transform the variant data into BED format
var.data = var.data %>% mutate(start = pos - 1, stop = pos) %>% select(chr, start, stop, Gene, Ind, Z, Outlier) %>%
	mutate(VarId = paste(chr, stop, Gene, sep = ':')) %>% bedr.sort.region(., verbose = F)

# Get list of promoter annotation files
prom.ann.files = list.files(paste0(dir, '/features/annotations/epigenomicsRoadmap/prom'), full.names = T)
er.tissues = gsub('\\.bed.*', '', gsub('.*/', '', prom.ann.files))
# For each promoter annotation BED file, intersect with the variation data
for(i in 1:length(prom.ann.files)){
	prom.bed.file = prom.ann.files[i]
	tissue = er.tissues[i]
	print(tissue)
	prom.bed = read.table(prom.bed.file, sep = '\t', header = F, stringsAsFactors = F)[, 1:3]
	colnames(prom.bed) = c('chr', 'start', 'stop')
	prom.bed = prom.bed %>% bedr.sort.region(., verbose = F)
	var.data = bedr(input = list(a = var.data, b = prom.bed), params = '-wa -c', method = 'intersect', verbose = F)
	names(var.data)[ncol(var.data)] = tissue
}
names(var.data) = c('chr', 'start', 'stop', 'Gene', 'Ind', 'Z', 'Outlier', 'VarId', er.tissues)
var.data[, c('Z', er.tissues)] = apply(var.data[, c('Z', er.tissues)], 2, function(x){as.numeric(as.character(x))})
var.data = var.data %>% mutate(chr = as.numeric(gsub('chr', '', chr)), Outlier = factor(Outlier, levels = out.states))
var.data$ERProms = rowSums(var.data[, er.tissues])

# Add allele frequency data from Xin
frq.data = read.table('/srv/scratch/restricted/goats/shared_data/fromXin/uk10k.gtex148sites.af.processed.txt', sep = '\t', header = F, stringsAsFactors = F)
colnames(frq.data) = c('chr', 'stop', 'NumAlleles', 'NumChroms', 'MajorAllele', 'MajorAlleleCount', 'MinorAlleles', 'MinorAlleleCount')
frq.data = frq.data %>% distinct(chr, stop, MajorAllele, MinorAlleles, .keep_all = T)

# Merge the UK10K allele count data with the rare variant data
var.data = merge(var.data, frq.data[, c(1:4, 6, 8)], all.x = T, by = c('chr', 'stop')) %>% mutate(VarId = paste(chr, stop, Gene, sep = ':'))

# Add annotation for over and under expression
var.data = var.data %>% mutate(Strat = ifelse(Z > 0, out.states.strat[3], out.states.strat[2])) %>% 
	mutate(Strat = factor(ifelse(Outlier == out.states.strat[1], out.states.strat[1], Strat), levels = out.states.strat))

# Reduce variants with multiple occurrences to a single row
var.data = var.data %>% select(-which(colnames(var.data) %in% c('Ind', 'Z'))) %>% distinct(VarId, Outlier, .keep_all = T)

# Stratify by gene so that over- and uner-expression rare variants will have gene matched controls
# Remove variants for genes without outliers
gene.states = var.data %>% filter(Outlier == out.states[2]) %>% select(Gene, Strat) %>% distinct(Gene, .keep_all = T)
colnames(gene.states)[2] = 'GeneStrat'
var.data = merge(var.data, gene.states, all.x = T, by = 'Gene') %>% mutate(GeneStrat = factor(GeneStrat, levels = out.states.strat[2:3])) %>%
	filter(!is.na(GeneStrat))

# Add minor allele count and MAF data
num.chroms = max(var.data$NumChroms, na.rm = T)
var.data = var.data %>% mutate(MAC = ifelse(is.na(MajorAlleleCount), num.chroms, MajorAlleleCount),
	mAC = ifelse(is.na(MinorAlleleCount), 0, MinorAlleleCount), 
	MAF = (mAC + 1) / (mAC + MAC + 2))

# Restrict to outliers with variants in promoter regions
prom.var.data = var.data %>% filter(ERProms > 0) %>% group_by(Gene) %>% filter(length(unique(Outlier)) > 1)

# Make the main figure plot: Figure 4A
prom.var.data = prom.var.data %>% mutate(mACPretty = ifelse(mAC > 4, '>4', paste0(mAC)))
prom.var.table = table(prom.var.data$mACPretty, prom.var.data$Outlier)
prom.var.table = prom.var.table %*% diag(1 / colSums(prom.var.table))
colnames(prom.var.table) = out.states
prom.var.df = as.data.frame(prom.var.table) %>% mutate(mAC = factor(rownames(.), levels = c(paste0(0:4), '>4'))) %>%
	melt(id.vars = 'mAC', variable.name = 'Outlier', value.name = 'Freq')

uk10k.prom.wleg = ggplot(data = prom.var.df, aes(x = mAC, y = Freq, fill = Outlier)) +
	geom_bar(stat = 'identity', position = 'dodge', colour = 'black') + theme_bw() + xlab('Minor allele count in UK10K') +
	ylab('Proportion of promoter variants') + scale_fill_manual(values = out.colors) + theme(legend.title=element_blank()) +
	theme(legend.key = element_blank())
uk10k.prom.wleg

# Use Wilcoxon rank sum test to obtain significance level
uk10k.sig.test = wilcox.test(MAF ~ Outlier, data = prom.var.data)
uk10k.sig.test
# W = 50550, p-value = 0.005974

# Save workspace image
save.image(paste0(dir, '/data/figure4a.uk10k.RData'))


