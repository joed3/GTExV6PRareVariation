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

# Function to aggregate replication Z-scores
agg.zs <- function(x){
	return(median(unique(sort(x))))
}

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

# Read in variant prioritization file from Yung
vars = read.table(paste0(dir, '/data/CRISPR/candidate_SNVs_120116.txt'), sep = '\t', header = T, stringsAsFactors = F)
vars = vars %>% mutate(Outlier = paste(gene, indiv, sep = ':'), sv = ifelse(is.na(sv), 0, sv), indel = ifelse(is.na(indel), 0, indel), 
		InGeuvadis = ifelse(is.na(InGeuvadis), 0, InGeuvadis), InSardinia = ifelse(is.na(InSardinia), 0, InSardinia))
dim(vars)
# 69157 (Gene, Individual, Variant) tuples

# Read in set of multi-tissue outliers
outliers = read.table(paste0(dir, '/data/outliers_medz_picked.txt'), sep = '\t', header = T, stringsAsFactors = F)
outliers = outliers %>% mutate(Outlier = paste(GENE, INDS, sep = ':'))

# Filter for annotated variants occurring at multi-tissue outliers
vars = filter(vars, Outlier %in% outliers$Outlier)
dim(vars)
# 734 (Gene, Individual, Variant) tuples

# Filter for variants without nearby SVs or indels
vars = filter(vars, !(sv | indel))
dim(vars)
# 553 (Gene, Individual, Variant) tuples

# Filter for coding variants
vars$Coding = apply(vars[, coding.annotations], 1, any)
vars = filter(vars, Coding)
dim(vars)
# 116 (Gene, Individual, Variant) tuples

# Add K562 expression data from ENCODE
k562.cols = c('gene_id', 'FPKM')

k562.cshl1 = read.table(paste0(dir, '/data/K562/ENCFF104VTJ_CSHL_1.tsv'), sep = '\t', header = T, stringsAsFactors = F)[, k562.cols]
colnames(k562.cshl1)[2] = 'FPKM.CSHL1'

k562.cshl2 = read.table(paste0(dir, '/data/K562/ENCFF201HGA_CSHL_2.tsv'), sep = '\t', header = T, stringsAsFactors = F)[, k562.cols]
colnames(k562.cshl2)[2] = 'FPKM.CSHL2'

k562.uconn1 = read.table(paste0(dir, '/data/K562/ENCFF553DDU_UConn_1.tsv'), sep = '\t', header = T, stringsAsFactors = F)[, k562.cols]
colnames(k562.uconn1)[2] = 'FPKM.UConn1'

k562.uconn2 = read.table(paste0(dir, '/data/K562/ENCFF811VBA_UConn_2.tsv'), sep = '\t', header = T, stringsAsFactors = F)[, k562.cols]
colnames(k562.uconn2)[2] = 'FPKM.UConn2'

k562 = merge(k562.cshl1, k562.cshl2, all = T, by = k562.cols[1])
k562 = merge(k562, k562.uconn1, all = T, by = k562.cols[1])
k562 = merge(k562, k562.uconn2, all = T, by = k562.cols[1])

k562$FPKM.Mean = rowMeans(k562[, -1], na.rm = T)

rownames(k562) = gsub('\\.[0-9]*', '', k562$gene_id)

vars[, colnames(k562)[-1]] = k562[gsub('\\.[0-9]*', '', vars$gene), -1]

# Summarize coding variant types
# How many of each type?
coding.counts = unlist(lapply(1:length(coding.annotations), function(x){sum(vars[, coding.annotations[x]])}))
coding.counts = data.frame(Annotations = coding.annotations.pretty, Count = coding.counts, stringsAsFactors = F)
coding.counts = coding.counts[order(coding.counts$Count), ]
coding.counts$Annotations = factor(coding.counts$Annotations, levels = coding.counts$Annotations)

coding.bar = ggplot(data = coding.counts, aes(x = Annotations, y = Count)) + geom_bar(aes(fill = Annotations), stat = 'identity', alpha = .6) + 
	theme_bw() + fontsizes + xlab('') + ylab('Count') + coord_flip() + scale_fill_manual(values = annot.colors) + guides(fill = F) + 
	annotate('text', x = coding.counts$Annotations, y = coding.counts$Count + 5, label = coding.counts$Count) +
	ggtitle('Coding annotations')

# How many types for each variant?
var.counts = data.frame(Count = rowSums(vars[, coding.annotations]))
var.counts.hist = as.vector(table(var.counts$Count))
var.counts.bar = ggplot(data = var.counts, aes(x = Count)) + geom_histogram(fill = 'dodgerblue4', alpha = .6) + theme_bw() + fontsizes +
	xlab('Number of coding annotations') + ylab('Number of variants') + 
	annotate('text', x = 1:max(var.counts$Count), y = var.counts.hist + 5, label = var.counts.hist) +
	ggtitle('Annotation overlap')

# Plot distribution of RIVER scores on gene and nucleotide level
river.nuc = ggplot(data = vars, aes(x = pFRgivenGnE)) + geom_histogram(fill = 'dodgerblue4', alpha = .6) + theme_bw() +
	fontsizes + xlab('RIVER score (nucleotide-level)') + ylab('Count') + geom_vline(xintercept = median(vars$pFRgivenGnE)) +
	annotate('text', x = .5, y = 20, label = paste('Median = ', round(median(vars$pFRgivenGnE), 3), sep = '')) +
	ggtitle('RIVER')

# Plot distribution of number of tissues
ntiss = ggplot(data = vars, aes(x = NObsTissue)) + geom_histogram(fill = 'dodgerblue4', alpha = .6) + theme_bw() + fontsizes +
	ylab('Count') + xlab('Number of tissues') + geom_vline(xintercept = median(vars$NObsTissue)) +
	annotate('text', x = 22, y = 10, label = paste('Median = ', round(median(vars$NObsTissue), 3), sep = '')) +
	ggtitle('Number of tissues')

# Plot relationship between number of tissues and effect size
ntiss.effect.corr = cor(vars$NObsTissue, abs(vars$medianZ), method = 'p', use = 'complete')

ntiss.effect = ggplot(data = vars, aes(x = NObsTissue, y = abs(medianZ))) + 
	geom_point(colour = 'dodgerblue4', alpha = .6, size = 3) + theme_bw() + 
	fontsizes + ylab('|Median Z-score|') + xlab('Number of tissues') + 
	annotate('text', x = 20, y = 9, label = 'Pearson~R^{2}==0.002', parse = T) +
	ggtitle('Tissues v. Effect size')

# Plot distribution of number of cis-eQTLs
neqtls = ggplot(data = vars, aes(x = NTopeQTLs)) + geom_histogram(fill = 'dodgerblue4', alpha = .6) + theme_bw() + fontsizes + 
	xlab('Number of cis-eQTLs') + ylab('Count') + geom_vline(xintercept = median(vars$NTopeQTLs)) +
	annotate('text', x = 10, y = 40, label = paste('Median = ', round(median(vars$NTopeQTLs), 3), sep = '')) +
	ggtitle('Number of cis-eQTLs')

# Plot relationship between effect size and number of cis-eQTLs
neqtls.effect.corr = cor(vars$NTopeQTLs, abs(vars$medianZ), use = 'complete', method = 'p')
neqtls.effect = ggplot(data = vars, aes(x = NTopeQTLs, y = abs(medianZ))) + geom_point(colour = 'dodgerblue4', alpha = .6, size = 3) +
	theme_bw() + fontsizes + xlab('Number of cis-eQTLs') + ylab('|Median Z-score|') + 
	annotate('text', x = 10, y = 8, label = 'Pearson~R^{2}==-0.098', parse = T) +
	ggtitle('cis-eQTLs v. Effect size')

# Make plot of mean FPKM level in K562 cells
k562.hist = ggplot(data = vars, aes(x = log10(FPKM.Mean + 1))) + geom_histogram(fill = 'dodgerblue4', alpha = .6) + theme_bw() + fontsizes +
	xlab('log10(K562 FPKM + 1)') + ylab('Count') + geom_vline(xintercept = median(log10(vars$FPKM.Mean + 1), na.rm = T)) +
	annotate('text', x = 2, y = 10, label = paste('Median = ', median(vars$FPKM.Mean, na.rm = T))) +
	ggtitle('K562 expression')

# Make figure with these summary plots
grid.arrange(coding.bar, 
	var.counts.bar,
	k562.hist, 
	ntiss, 
	ntiss.effect, 
	neqtls, 
	neqtls.effect,
	river.nuc, 
	ncol = 3, nrow = 3)

# Order the variants by nucleotide-level RIVER score
vars.ordered = vars[order(vars$pFRgivenGnE, decreasing = T), ]

# Write out list of variants
write.table(vars.ordered[, -c(35:37)], paste0(dir, '/data/CRISPR/prioritized.variants.for.crispr.txt'), sep = '\t', col.names = T, row.names = F, quote = F)




