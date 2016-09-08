#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(ggplot2)
require(reshape2)
require(plyr)

#------------- FUNCTIONS

# Function to convert Xin's colors to RGB strings
rgb2string <- function(x){
	return(rgb(red = x['R'], green = x['G'], blue = x['B']))
}

# Aggregation function for unmelting Xin's data
# Returns 0 or 1 based on whether the annotation exists for the individual or not 
# (i.e. whether the individual has a rare variant for the gene)
cat.agg <- function(x){
	return((length(x) > 0) + 0)
}

# Function to compute log odds ratios and CIs for Xin's annotation
# Return DF of estimates ordered by p-value with SV annotations appearing first
vep.estimate <- function(data, features, sv.annots, scale = T){
	estims = c()
	for(i in 1:length(features)){
		if(scale){
			mod = summary(glm(isOutlier ~ scale(data[, features[i]]), data = data, family = binomial))$coefficients
			estims.temp = data.frame(Cat = features[i], Estim = mod[2, 1], Std = mod[2, 2], Pval = mod[2, 4], stringsAsFactors = F)	
			estims = rbind(estims, estims.temp)
		}
		else{
			mod = summary(glm(isOutlier ~ data[, features[i]], data = data, family = binomial))$coefficients
			estims.temp = data.frame(Cat = features[i], Estim = mod[2, 1], Std = mod[2, 2], Pval = mod[2, 4], stringsAsFactors = F)	
			estims = rbind(estims, estims.temp)
		}
	}
	sv.estims = estims[estims$Cat %in% sv.annots, ]
	sv.estims = sv.estims[order(sv.estims$Pval), ]
	other.estims = estims[!(estims$Cat %in% sv.annots), ]
	other.estims = other.estims[order(other.estims$Pval), ]
	estims = rbind(sv.estims, other.estims)
	estims$Cat = factor(as.character(estims$Cat), levels = as.character(estims$Cat))
	return(estims)
}

#------------- MAIN

# Set median Z-score outlier threshold
out.thresh = 2

# Read in file with Xin's colors
colors = read.csv(paste0(dir, '/data/fromXin/colors.from.xin.txt'), sep = '\t', header = T, stringsAsFactors = F)

# Define features to test
comp.levels = colors$HumanID

# Make final color vector for plotting
colors.final = colors$hex
names(colors.final) = colors$HumanID

# Define SV annotations
sv.annots = names(colors.final)[1:5]

# Read in Xin's file of outliers with rare variant categories
# Map rare variant categories to human readable format
outlier.cats.raw = read.table(paste0(dir, '/data/fromXin/Fig3a_outlierVariants_filtered.txt'), header = T, sep = '\t')
outlier.cats.raw$compoundId = mapvalues(outlier.cats.raw$compoundId, colors$compoundId, colors$HumanID)

# Unmelt the data to have a single column for each categories
outlier.cats = dcast(data = outlier.cats.raw, GeneId + INDid + isOutlier ~ compoundId, fun.aggregate = cat.agg)

# Calculate enrichments for each categories
# First with scaling and then without
estims.scaled = vep.estimate(outlier.cats, comp.levels, sv.annots, scale = T)

# Make forest plot (vertical)
estims.scaled$Cat = factor(as.character(estims.scaled$Cat), levels = rev(levels(estims.scaled$Cat)))

outlier.cats.forest.scaled = ggplot(data = estims.scaled[!(estims.scaled$Cat == 'No variant'), ], aes(x = Cat, y = Estim)) +
    geom_pointrange(aes(x = Cat, ymin = Estim - 1.96 * Std, ymax = Estim + 1.96 * Std, colour = Cat)) +
    geom_vline(xintercept = 7.5, linetype = "dashed", colour = "darkgrey") +
    geom_hline(yintercept = 0) +
    theme_bw() + xlab('') + ylab('Log odds ratio') + scale_colour_manual(values = colors.final) + guides(colour = F) +
    coord_flip() + ylim(c(-0.05, .3))
outlier.cats.forest.scaled

# Save workspace image
save.image(paste0(dir, '/data/figure3a.rare.variant.class.enrichments.RData'))



