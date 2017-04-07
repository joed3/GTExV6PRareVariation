#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(ggplot2)
require(reshape2)
require(data.table)
require(scales)
require(plyr)

#------------- FUNCTIONS

source('enrichment.functions.R')

# Function to impute missing TFBS features from CADD to 0
impute.tfbs <- function(features){
	features$any_tfbs_CADD = ifelse(is.na(features$any_tfbs_CADD), 0, features$any_tfbs_CADD)
	features$any_tfbs_CADD_peaks = ifelse(is.na(features$any_tfbs_CADD_peaks), 0, features$any_tfbs_CADD_peaks)
	features$max_tfbs_CADD = ifelse(is.na(features$max_tfbs_CADD), 0, features$max_tfbs_CADD)
	features$max_tfbs_CADD_peaks = ifelse(is.na(features$max_tfbs_CADD_peaks), 0, features$max_tfbs_CADD_peaks)
	features$max_tfbs_CADD_peaks_max = ifelse(is.na(features$max_tfbs_CADD_peaks_max), 0, features$max_tfbs_CADD_peaks_max)
	return(features)
}

#------------- MAIN

# Define estimate types
types = c('10kb with PC', '200kb with PC')

# Read in features for indels and SNVs
# 10kb distance threshold with PC
features.snvs.10kb = read.table(paste0(dir, '/data/medz/10kb_genebody_features_SNPs_MAF0-1.txt'), sep = '\t', header = T, stringsAsFactors = F)

# 200kb distance threshold with PC
features.snvs.200kb = read.table(paste0(dir, '/data/medz/200kb_genebody_features_SNPs_MAF0-1.txt'), sep = '\t', header = T, stringsAsFactors = F)

# Replace missing TFBS features from CADD with 0
features.10kb = impute.tfbs(features.snvs.10kb)
features.200kb = impute.tfbs(features.snvs.200kb)

# Find the estimated log odds estimates and CIs
# Scaled
estims.10kb.scaled = single_logit_core(features.10kb)
estims.10kb.scaled$TYPE = factor(types[1], levels = types)

estims.200kb.scaled = single_logit_core(features.200kb)
estims.200kb.scaled$TYPE = factor(types[2], levels = types)

# Replace enhancer LORs in 10kb with PC with those from 200kb with PC
estims.10kb.scaled[estims.10kb.scaled$FEAT == 'any_enhancer', ] = estims.200kb.scaled[estims.200kb.scaled$FEAT == 'any_enhancer', ]
estims.10kb.scaled[estims.10kb.scaled$FEAT == 'max_enhancer', ] = estims.200kb.scaled[estims.200kb.scaled$FEAT == 'max_enhancer', ]

# Remove dyadic features and n_variants since they have no variation
dynames = c('any_dyadic', 'max_dyadic', 'n_variants')

estims.10kb.scaled = estims.10kb.scaled[!(estims.10kb.scaled$FEAT %in% dynames), ]

# Make plot of additive effects
# Organize feature list into similar groups and remove redundant annotations (i.e. choose any over max when possible)
add.names.raw = c('closest_variant', 
	'any_promoter', 
	'any_enhancer', 
	'any_tfbs_CADD', 
	'max_CpG', 
	'max_cadd', 
	'max_ver_phylop',
	'max_ver_phastCons', 
	'max_fitCons', 
	'max_GerpN', 
	'max_GerpS')

add.names = c('Distance to TSS', 
	'Promoter', 
	'Enhancer', 
	'TFBS', 
	'CpG', 
	'CADD', 
	'PhyloP', 
	'PhastCons', 
	'fitCons', 
	'GerpN', 
	'GerpS')

estims.10kb.scaled = estims.10kb.scaled[estims.10kb.scaled$FEAT %in% add.names.raw, ]
estims.10kb.scaled$FEAT = factor(as.character(estims.10kb.scaled$FEAT), levels = add.names.raw)
estims.10kb.scaled$FEAT = mapvalues(estims.10kb.scaled$FEAT, from = add.names.raw, to = add.names)

estims.10kb.scaled = estims.10kb.scaled[nrow(estims.10kb.scaled):1, ]
estims.10kb.scaled$FEAT = factor(as.character(estims.10kb.scaled$FEAT), levels = add.names[length(add.names):1])

add.colors = c('mediumorchid3', rep('mediumpurple4', 4), 'tomato', rep('maroon3', 5))
add.colors = add.colors[length(add.colors):1]

add.effects.plot.scaled = ggplot(data = estims.10kb.scaled, aes(x = FEAT, y = BETA, colour = FEAT)) + 
	geom_pointrange(aes(x = FEAT, ymin = BETA - 1.96 * STD, ymax = BETA + 1.96 * STD)) +
	theme_bw() + xlab('') + ylab('Log odds ratio') + 
        geom_vline(xintercept = c(5.5, 6.5, 10.5), linetype = "dashed", colour = "darkgrey") +
        geom_hline(yintercept = 0) +
	scale_colour_manual(values = add.colors) + guides(colour = F) +
	theme(strip.text.x = element_text(face='bold')) + coord_flip()
add.effects.plot.scaled

# Save workspace image
save(add.effects.plot.scaled, file = paste0(dir, '/data/figure3b.feature.enrichments.RData'))
