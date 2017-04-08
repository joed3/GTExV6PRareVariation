#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(ggplot2)
require(reshape2)
require(data.table)
require(scales)
require(dplyr)
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
window.types = c('10kb with PC', '200kb with PC')
outlier.types = c('Underexpression', 'Overexpression')

# Read in features for indels and SNVs
# 10kb distance threshold with PC
features.snvs.10kb = read.table(paste0(dir, '/data/medz/10kb_genebody_features_SNPs_MAF0-1.txt'), sep = '\t', header = T, stringsAsFactors = F)
# 200kb distance threshold with PC
features.snvs.200kb = read.table(paste0(dir, '/data/medz/200kb_genebody_features_SNPs_MAF0-1.txt'), sep = '\t', header = T, stringsAsFactors = F)

# Remove dyadic features and n_variants since they have no variation
# Only focus on the enhancer features for the 200kb df
dynames = c('any_dyadic', 'max_dyadic', 'n_variants')
enhnames = c('any_enhancer', 'max_enhancer')

# Replace missing TFBS features from CADD with 0
features.10kb = impute.tfbs(features.snvs.10kb) %>% select(which(!(colnames(features.snvs.10kb) %in% dynames)))
features.200kb = impute.tfbs(features.snvs.200kb) %>% select(1, 29:31, which(colnames(features.snvs.200kb) %in% enhnames))

# Split over and under expression outliers
f10.over = subset(features.10kb, RANK > 0)
f10.under = subset(features.10kb, RANK < 0)
f200.over = subset(features.200kb, RANK > 0)
f200.under = subset(features.200kb, RANK < 0)

# Find the estimated log odds estimates and CIs
# Scaled
estims.10kb.scaled = single_logit_core(features.10kb) %>% mutate(TYPE = factor(window.types[1], levels = window.types))
estims.200kb.scaled = single_logit_core(features.200kb) %>% mutate(TYPE = factor(window.types[2], levels = window.types))

estims.10kb.over = single_logit_core(f10.over) %>% mutate(TYPE = factor(outlier.types[2], levels = outlier.types))
estims.10kb.under = single_logit_core(f10.under) %>% mutate(TYPE = factor(outlier.types[1], levels = outlier.types))
estims.200kb.over = single_logit_core(f200.over) %>% mutate(TYPE = factor(outlier.types[2], levels = outlier.types))
estims.200kb.under = single_logit_core(f200.under) %>% mutate(TYPE = factor(outlier.types[1], levels = outlier.types))

# Replace enhancer LORs in 10kb with PC with those from 200kb with PC
estims.10kb.scaled[estims.10kb.scaled$FEAT == 'any_enhancer', ] = estims.200kb.scaled[estims.200kb.scaled$FEAT == 'any_enhancer', ]

estims.10kb.over[estims.10kb.over$FEAT == 'any_enhancer', ] = estims.200kb.over[estims.200kb.over$FEAT == 'any_enhancer', ]
estims.10kb.under[estims.10kb.under$FEAT == 'any_enhancer', ] = estims.200kb.under[estims.200kb.under$FEAT == 'any_enhancer', ]


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

estims.10kb.scaled = filter(estims.10kb.scaled, FEAT %in% add.names.raw) %>%
    mutate(FEAT = factor(mapvalues(FEAT, from = add.names.raw, to = add.names), levels = add.names[length(add.names):1]))

estims.10kb.over = filter(estims.10kb.over, FEAT %in% add.names.raw) %>%
    mutate(FEAT = factor(mapvalues(FEAT, from = add.names.raw, to = add.names), levels = add.names[length(add.names):1]))
estims.10kb.under = filter(estims.10kb.under, FEAT %in% add.names.raw) %>%
            mutate(FEAT = factor(mapvalues(FEAT, from = add.names.raw, to = add.names), levels = add.names[length(add.names):1]))
estims.10kb = rbind(estims.10kb.over, estims.10kb.under)

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

type.shapes = c(25, 24)
names(type.shapes) = outlier.types

add.effects.plot.over.under = ggplot(data = estims.10kb, aes(x = FEAT, y = BETA, colour = FEAT, shape = TYPE)) +
    geom_pointrange(aes(x = FEAT, ymin = BETA - 1.96 * STD, ymax = BETA + 1.96 * STD), position = position_dodge(width = .5)) +
    theme_bw() + xlab('') + ylab('Log odds ratio') +
    geom_vline(xintercept = c(5.5, 6.5, 10.5), linetype = "dashed", colour = "darkgrey") +
    geom_hline(yintercept = 0) +
    scale_colour_manual(values = add.colors) + guides(colour = F) +
    theme(strip.text.x = element_text(face='bold')) + coord_flip() +
    scale_shape_manual(values = type.shapes, name = '', guide = guide_legend(reverse = T))
add.effects.plot.over.under

# Save workspace image
save(add.effects.plot.scaled, add.effects.plot.over.under, file = paste0(dir, '/data/figure3b.feature.enrichments.RData'))
