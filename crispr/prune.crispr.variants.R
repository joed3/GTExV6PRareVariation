#!/usr/bin/env Rscript

rm(list = ls())

dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(data.table)
require(reshape2)
require(plyr)
require(dplyr)
require(cowplot)
require(data.table)
require(gridExtra)

#-------------- FUNCTIONS


#-------------- MAIN

# Coding variant annotations and colors
coding.annotations = c("StopGained", 
	"StopLost",                    
	"SpliceAcceptorVariant",
	"SpliceDonorVariant",        
	"StartLost",
	"MissenseVariant",             
	"SpliceRegionVariant",
	"StopRetainedVariant",           
	"SynonymousVariant",
	"CodingSequenceVariant",         
	"X5PrimeUTRVariant",   
	"X3PrimeUTRVariant")

coding.annotations.pretty = c("Stop Gained", 
	"Stop Lost",                    
	"Splice Acceptor Variant",
	"Splice Donor Variant",        
	"Start Lost",
	"Missense Variant",             
	"Splice Region Variant",
	"Stop Retained Variant",           
	"Synonymous Variant",
	"Coding Sequence Variant",         
	"5' Prime UTR Variant",   
	"3' Prime UTR Variant")

annot.colors = c('#053061', 
	'#053061', 
	'#2166AC', 
	'#2166AC', 
	'#053061', 
	'#FDE595', 
	'#2166AC', 
	'#FDE595', 
	'#FDE595', 
	'#FDE595', 
	'#4393C3', 
	'#4393C3')
names(annot.colors) = coding.annotations.pretty

noncoding.annotations = c("NonCodingTranscriptExonVariant", 
	"NMDTranscriptVariant", 
	"NonCodingTranscriptVariant", 
	"IntronVariant")

noncoding.annotations.pretty = c("Non-Coding Transcript Exon Variant", 
	"NMD Transcript Variant", 
	"Non-Coding Transcript Variant", 
	"Intron Variant")

# Plotting constants
fontsize = 11
fontsizes = theme(axis.text = element_text(size = fontsize),
                  axis.title = element_text(size = fontsize),
                  strip.text = element_text(size = fontsize + 1),
                  legend.text = element_text(size = fontsize - 1),
                  legend.title = element_text(size = fontsize - 1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  legend.key = element_blank())

# Read in RIVER data
river = read.table(paste0(dir, '/data/RIVER.Posteriors.All.Final.V6P.111416.txt'), sep = '\t', header = T, stringsAsFactors = F)

# Calculate RIVER percentiles
# Use the .995 quantile as threshold
river.quants = quantile(river$pZgivenGE, prob = seq(0, 1, .001))
river.thresh = river.quants[length(river.quants) - 5]

# Read in the list of prioritzed variants
vars = read.table(paste0(dir, '/data/CRISPR/prioritized.variants.for.crispr.txt'), sep = '\t', header = T, stringsAsFactors = F)
dim(vars)
# 116 variants

# Select variants that are in the top 99.9% RIVER score, have mean FPKM > 10, and have effect size > 4
effect.thresh = 4
fpkm.thresh = 10

# Filter on expression level in K562
vars.pared = filter(vars, (!(is.na(FPKM.Mean)) & FPKM.Mean > fpkm.thresh))
dim(vars.pared)
# 59 variants

# Filter on effect size
vars.pared = filter(vars.pared, absmedZ > effect.thresh)
vars.pared[vars.pared$absmedZ > effect.thresh, ]
dim(vars.pared)
# 13 variants

# Filter on RIVER score
vars.pared = vars.pared[vars.pared$pFRgivenGgE > river.thresh, ]
dim(vars.pared)
# 13 variants

# Summarize coding variant types
# How many of each type?
coding.counts = unlist(lapply(1:length(coding.annotations), function(x){sum(vars.pared[, coding.annotations[x]])}))
coding.counts = data.frame(Annotations = coding.annotations.pretty, Count = coding.counts, stringsAsFactors = F)
coding.counts = coding.counts[order(coding.counts$Count), ]
coding.counts$Annotations = factor(coding.counts$Annotations, levels = coding.counts$Annotations)

coding.bar = ggplot(data = coding.counts, aes(x = Annotations, y = Count)) + geom_bar(aes(fill = Annotations), stat = 'identity', alpha = .6) + 
	theme_bw() + fontsizes + xlab('') + ylab('Count') + coord_flip() + scale_fill_manual(values = annot.colors) + guides(fill = F) + 
	annotate('text', x = coding.counts$Annotations, y = coding.counts$Count + .1, label = coding.counts$Count) +
	ggtitle('Coding annotations')

# How many types for each variant?
var.counts = data.frame(Count = rowSums(vars.pared[, coding.annotations]))
var.counts.hist = as.vector(table(var.counts$Count))
var.counts.bar = ggplot(data = var.counts, aes(x = Count)) + geom_histogram(fill = 'dodgerblue4', alpha = .6) + theme_bw() + fontsizes +
	xlab('Number of coding annotations') + ylab('Number of variants') + 
	annotate('text', x = 1:max(var.counts$Count), y = var.counts.hist + .5, label = var.counts.hist) +
	ggtitle('Annotation overlap')

grid.arrange(coding.bar, var.counts.bar, ncol = 2, nrow = 1)

# Simplify consequences from Yung into a single column
vars.pared$CodingCSQ = gsub("'", '', gsub(' ', '_', unlist(apply(vars.pared[, coding.annotations], 1, function(x){paste(coding.annotations.pretty[which(x == 1)], collapse = ';')}))))
vars.pared$CodingCSQ = ifelse(vars.pared$CodingCSQ == '', NA, vars.pared$CodingCSQ)
vars.pared$NoncodingCSQ = gsub("'", '', gsub(' ', '_', unlist(apply(vars.pared[, noncoding.annotations], 1, function(x){paste(noncoding.annotations.pretty[which(x == 1)], collapse = ';')}))))
vars.pared$NoncodingCSQ = ifelse(vars.pared$NoncodingCSQ == '', NA, vars.pared$NoncodingCSQ)
vars.pared = vars.pared[, -c(19:34)]

# Add ref/alt allele and genotype info
# Add CADD scores and 1KG MAFs
# Reformat to VCF
inds = unique(sort(vars.pared$indiv))
vars.pared$VariantId = paste(vars.pared$chr, vars.pared$pos, sep = ':')
vars.pared[, c('MAF', 'Genotype', 'REF', 'ALT', 'CADD.Raw', 'CADD.Phred', 'EAS', 'AMR', 'AFR', 'EUR', 'SAS')] = NA
for(i in 1:length(inds)){
	ind = inds[i]
	snps = read.table(paste0(dir, '/features/variantBeds/individuals/', ind, '_SNPs.bed.gz'), sep = '\t', header = F, stringsAsFactors = F)[, -2]
	colnames(snps) = c('chr', 'pos', 'MAF', 'Genotype', 'Ref', 'Alt')
	snps$chr = as.numeric(gsub('chr', '', snps$chr))
	snps$VariantId = paste(snps$chr, snps$pos, sep = ':')
	rownames(snps) = snps$VariantId
	snp.features = read.table(paste0(dir, '/features/bySite/', ind, '_SNPs_features.bed.gz'), sep = '\t', header = T, stringsAsFactors = F)
	snp.features[, 1] = gsub('chr', '',snp.features[, 1])
	rownames(snp.features) = paste(snp.features[, 1], snp.features[, 3], sep = ':')
	indices = which(vars.pared$indiv == ind)
	for(j in 1:length(indices)){
		vars.pared[indices[j], c('MAF', 'Genotype', 'REF', 'ALT')] = snps[vars.pared$VariantId[indices[j]], c('MAF', 'Genotype', 'Ref', 'Alt')]
		vars.pared[indices[j], c('CADD.Raw', 'CADD.Phred', 'EAS', 'AMR', 'AFR', 'EUR', 'SAS')] = snp.features[vars.pared$VariantId[indices[j]], c(6:7, 27:31)]
	}
}
vars.pared[, c('QUAL', 'FILTER', 'INFO')] = NA
vars.pared = vars.pared[, c('chr', 
	'pos', 
	'VariantId', 
	'REF', 
	'ALT', 
	'QUAL', 
	'FILTER', 
	'INFO', 
	'Genotype', 
	'MAF', 
	'indiv', 
	'gene', 
	'medianZ', 
	'CodingCSQ', 
	'NoncodingCSQ', 
	colnames(vars.pared)[c(5:11, 14:15, 18:25, 33:39)])]
colnames(vars.pared)[1] = '#CHROM'
colnames(vars.pared)[2] = 'POS'
colnames(vars.pared)[3] = 'ID'

# Add MAFs from UK10K
uk10k = read.table('/users/xli6/data/xin/gtex/replication/gtex_allrare.txt', sep = '\t', header = T, stringsAsFactors = F)[, c(1:2, 6:7, 15)]
uk10k$INDid = paste0('GTEX-', uk10k$INDid)
uk10k$ID = paste(uk10k$variantId_1, uk10k$variantId_2, uk10k$GeneId, uk10k$INDid, sep = ':')
uk10k = uk10k[order(uk10k$ID), ]
uk10k = unique(uk10k)
rownames(uk10k) = uk10k$ID

N.UK10K = 7562
vars.pared$UK10K = uk10k[paste(vars.pared$ID, vars.pared$gene, vars.pared$indiv, sep = ':'), 'ALTc_uk10k_7562'] / (2 * N.UK10K)

#---- Add expression and Z-score info to the CRISPR variants
genes = unique(sort(vars.pared$gene))

# Read in RPKM flat file 
# Set Gene as the key
rpkm = fread(paste0(dir, '/preprocessing/gtex_2015-01-12_rpkm.txt'))
setkey(rpkm, Gene)

# Subset RPKM data to the genes we are considering for CRISPR
rpkm = rpkm[Gene %in% genes]

# Read in Z-scores file
# Set Gene as key
scores = fread(paste0(dir, '/preprocessing/gtex_2015-01-12_normalized_expression.txt'))
setkey(scores, Gene)

# Subset Z-scores for the genes we are considering for CRISPR
scores = scores[Gene %in% genes]

# Define the set of tissues available for analysis
tissues = unique(sort(scores$Tissue))

# Make plots showing RPKM and Z-scores for each (Gene, Individual, Variant) tuple
rpkm.outliers = matrix(, ncol = length(tissues), nrow = nrow(vars.pared))
colnames(rpkm.outliers) = tissues

fc.outliers = matrix(, ncol = length(tissues), nrow = nrow(vars.pared))
colnames(fc.outliers) = tissues

zscore.outliers = matrix(, ncol = length(tissues), nrow = nrow(vars.pared))
colnames(zscore.outliers) = tissues

for(i in 1:nrow(vars.pared)){
	gene = vars.pared$gene[i]
	id = vars.pared$indiv[i]
	# RPKM
	rpkm.temp = rpkm[Gene == gene]
	rpkm.temp = rpkm.temp[which(!is.na(rpkm.temp[, id, with = F])), ]
	rpkm.temp.medians = apply(rpkm.temp[, 3:nrow(rpkm.temp), with = F], 1, median, na.rm = T)
	rpkm.temp.melted = melt(rpkm.temp, id.vars = c('Tissue'), measure.vars = names(rpkm.temp)[3:ncol(rpkm.temp)], value.name = 'RPKM', variable.name = 'INDS')
	rpkm.temp.melted$Tissue = factor(rpkm.temp.melted$Tissue, levels = rev(unique(sort(rpkm.temp.melted$Tissue))))
	rpkm.outlier = rpkm.temp.melted[rpkm.temp.melted$INDS == id, ]
	rpkm.outliers[i, as.character(rpkm.outlier$Tissue)] = rpkm.outlier$RPKM
	fc.outliers[i, as.character(rpkm.outlier$Tissue)] = rpkm.outlier$RPKM / rpkm.temp.medians
	# Z-scores
	scores.temp = scores[Gene == gene]
	scores.temp = scores.temp[which(!is.na(scores.temp[, id, with = F])), ]
	scores.temp.melted = melt(scores.temp, id.vars = c('Tissue'), measure.vars = names(scores.temp)[3:ncol(scores.temp)], value.name = 'Zscore', variable.name = 'INDS')
	scores.temp.melted$Tissue = factor(scores.temp.melted$Tissue, levels = rev(unique(sort(scores.temp.melted$Tissue))))
	zscore.outlier = scores.temp.melted[scores.temp.melted$INDS == id, ]
	zscore.outliers[i, as.character(zscore.outlier$Tissue)] = zscore.outlier$Zscore
}

# Write out the RPKM and Z-score data for these variants
rpkm.outliers = cbind(vars.pared[, 1:12], as.data.frame(rpkm.outliers))
write.table(rpkm.outliers, paste0(dir, '/data/CRISPR/crispr.candidates.pruned.rpkm.vcf'), sep = '\t', col.names = T, row.names = F, quote = F)

fc.outliers = cbind(vars.pared[, 1:12], as.data.frame(fc.outliers))
write.table(fc.outliers, paste0(dir, '/data/CRISPR/crispr.candidates.pruned.fc.vcf'), sep = '\t', col.names = T, row.names = F, quote = F)

zscore.outliers = cbind(vars.pared[, 1:12], as.data.frame(zscore.outliers))
write.table(zscore.outliers, paste0(dir, '/data/CRISPR/crispr.candidates.pruned.zscores.vcf'), sep = '\t', col.names = T, row.names = F, quote = F)

# Write out this pruned set of variants with annotations
write.table(vars.pared, paste0(dir, '/data/CRISPR/crispr.candidates.pruned.vcf'), sep = '\t', col.names = T, row.names = F, quote = F)

# Generate BED file for these variants
radius = 2e5
vars.pared.bed = vars.pared[, c(1:5, 11:13, 20, 22)]
vars.pared.bed = cbind(vars.pared.bed[, 1], 
	vars.pared.bed[, 2] - 1 - radius, 
	vars.pared.bed[, 2] + radius, 
	vars.pared.bed[, 2], 
	vars.pared.bed[, 4:ncol(vars.pared.bed)])

write.table(vars.pared.bed[order(as.character(vars.pared.bed[, 1]), vars.pared.bed[, 2]), ], paste0(dir, '/data/CRISPR/crispr.candidates.pruned.bed'), sep = '\t', col.names = F, row.names = F, quote = F)





