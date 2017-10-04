#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(ggplot2)
require(reshape2)
require(plyr)

#-------------------- FUNCTIONS

# Function to convert Xin's colors to RGB strings
rgb2string <- function(x){
	return(rgb(red = x['R'], green = x['G'], blue = x['B']))
}

# Function to format annotation labels to fit on X-axis of effects plot
format.labels <- function(char.vec){
	for(i in 1:length(char.vec)){
		if(char.vec[i] == 'TSS non-coding'){
			char.vec[i] = 'TSS'
		}
		else if(char.vec[i] == 'Top 1% conserved non-coding'){
			char.vec[i] = 'Conserved'
		}
	}
	return(char.vec)
}

#-------------------- MAIN

# load necessary data
load(paste0(dir, '/data/figure3a.rare.variant.class.enrichments.RData'))

# Add row for non-outliers
estims.scaled = rbind(estims.scaled, data.frame(Cat = 'Non-outlier', Estim = NA, Std = NA, Pval = NA))

# Format the annotation labels to fit on the X-axis
estims.scaled$Cat = format.labels(as.character(estims.scaled$Cat))
colors$HumanID = format.labels(colors$HumanID)

# Read in Xin's effect size data
effects = read.table(paste0(dir, '/data/fromXin/Fig3d_effectsize_filtered.txt'), sep = '\t', header = T)
ase.effects = read.table(paste0(dir, '/data/fromXin/Fig3e_ASEeffectsize_filtered.txt'), sep = '\t', header = T)

# Define final colors vector
colors.final = c(colors$hex, 'darkgrey')
names(colors.final) = c(colors$HumanID, 'Non-outlier')

type.color = c('black', NA)
names(type.color) = c('Outlier', 'Non-outlier')

# Map levels of the compound ID to human readable form
effects$compoundId = mapvalues(effects$compoundId, from = c(colors$compoundId, 'nonOutlier'), to = c(colors$HumanID, 'Non-outlier'))
effects$compoundId = factor(as.character(effects$compoundId), levels = estims.scaled$Cat)
ase.effects$compoundId = mapvalues(ase.effects$compoundId, from = c(colors$compoundId, 'nonOutlier'), to = c(colors$HumanID, 'Non-outlier'))
ase.effects$compoundId = factor(as.character(ase.effects$compoundId), levels = estims.scaled$Cat)

# For each data frame, transform it to have the value of interest in a column (named value)
# and the name of that column in another column (names variable)
# then stack the two data frame, keeping the columns needed for plotting
effects.melted = melt(effects, id.vars = colnames(effects)[colnames(effects) != "medianZ"])
ase.effects.melted = melt(ase.effects, id.vars = colnames(ase.effects)[colnames(ase.effects) != "mean_ASEeffectSize"])
colsToKeep = c("compoundId", "variable", "value")
combined.effects = rbind(effects.melted[, colsToKeep], ase.effects.melted[, colsToKeep])

combined.effects$TYPE = factor(ifelse(combined.effects$compoundId == 'Non-outlier', 'Non-outlier', 'Outlier'), levels = names(type.color))

combined.effects$variable = factor(combined.effects$variable, labels = c("Median Z-score", "Mean ASE"))

median.nonoutlier.effect = median(effects.melted$value[effects.melted$compoundId == "Non-outlier"], na.rm = T)
median.nonoutlier.ase = median(ase.effects.melted$value[ase.effects.melted$compoundId == "Non-outlier"], na.rm = T)
median.data = data.frame(cat = c("Median Z-score", "Mean ASE across tissues"), yval = c(median.nonoutlier.effect, median.nonoutlier.ase), stringsAsFactors = F)
median.data$cat = factor(median.data$cat, levels = c("Median Z-score", "Mean ASE across tissues")) # match level order to combined.effects

# Plot effect sizes by annotation grouped by variant type
# For now, not including the other category. May want to include it in the supplement?
effects.plot = ggplot(data = combined.effects, aes(x = compoundId, y = value, group = compoundId)) +
    geom_hline(data = median.data, aes(yintercept = yval), colour = 'darkgrey') +
    geom_violin(aes(fill = compoundId), scale = 'width', width = .5) +
    geom_point(aes(colour = TYPE)) +
    geom_jitter(aes(colour = TYPE), width = .2) +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = 'left') +
    theme_bw() + guides(fill = F, colour = F) +
    scale_fill_manual(values = colors.final) +
    scale_colour_manual(values = type.color) +
    xlab('') + ylab('') +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(hjust = 0))
effects.plot

# Save workspace image
save(effects.plot, file = paste0(dir, '/data/figure3de.outlier.effect.size.RData'))
