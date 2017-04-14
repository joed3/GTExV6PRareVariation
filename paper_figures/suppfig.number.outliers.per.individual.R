#!/usr/bin/env Rscript

## PATHS TO SET ###############
dir = Sys.getenv('RAREVARDIR')
phenoFile = 'GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt' # adjust path as necessary

# Load required packages
require(ggplot2)
require(cowplot)


#---------------- MAIN

# Common axis font sizes
fontsize = 9
axisFontSizes = theme(axis.text = element_text(size = fontsize),
                      axis.title = element_text(size = fontsize),
                      legend.title = element_text(size = fontsize),
                      legend.text = element_text(size = fontsize),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())

# Read in phenotypes file
phenos = read.csv(phenoFile, sep = '\t', header = T, row.names = 1)

# Convert race and gender categories to human readable format
races = c('Asian', 'Black', 'White', 'Native\nAmerican')
races[98] = '    Unknown'
races[99] = races[98]
phenos$RACE = factor(races[phenos$RACE], levels = races[c(1:4, 98)])

genders = c('Male', 'Female')
phenos$GENDER = factor(genders[phenos$GENDER], levels = genders)

# Read in number of outliers per individual
ind.counts = read.table(paste0(dir, '/data/outliers_medz_picked_counts_per_ind.txt'), sep = '\t', header = F, stringsAsFactors = F)
colnames(ind.counts) = c('IND', 'COUNT')

# Sort individuals by number of outliers
ind.counts = ind.counts[order(ind.counts$COUNT), ]
ind.counts$IND = factor(ind.counts$IND, levels = ind.counts$IND)

# Define threshold for number of outliers per individual
count.thresh = 50

# Add column indicating individual was kept for further analyses or had too many outliers
ind.counts$Kept = as.factor(ifelse(ind.counts$COUNT < count.thresh, 'Good', 'Bad'))
kept.colors = c('darkgrey', 'dodgerblue3')
names(kept.colors) = c('Bad', 'Good')

# Make plot of the number of outliers per individual
ind.counts.plot = ggplot(data = ind.counts, aes(x = IND, y = COUNT)) + geom_bar(aes(fill = Kept, colour = Kept), stat = 'identity') + axisFontSizes + theme_bw() + 
	theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
	xlab(paste('Individuals (', nrow(ind.counts), ')', sep = '')) + ylab('Number of genes') + geom_abline(intercept = count.thresh, slope = 0) +
	scale_fill_manual(values = kept.colors) + scale_colour_manual(values = kept.colors) + guides(fill = F, colour = F) + 
        scale_y_continuous(breaks = c(0, 10, seq(20, max(ind.counts$COUNT), 30))) + theme(axis.text = element_text(size = fontsize),
                                                                                          axis.title = element_text(size = fontsize))
# Remove the outlier individuals
ind.counts = ind.counts[ind.counts$COUNT < count.thresh, ]

# Break the number of outliers down by race and gender
ind.counts[, c('RACE', 'GENDER', 'AGE', 'BMI', 'ISCH')] = phenos[ind.counts$IND, c('RACE', 'GENDER', 'AGE', 'BMI', 'TRISCHD')]

race.col = c('darkolivegreen3','darkgoldenrod2','dodgerblue2','#FF3D3D','darkgrey')

counts.by.race = ggplot(data = ind.counts, aes(x = RACE, y = COUNT)) +
    geom_violin(aes(fill = RACE), scale = 'width', width = .5) + theme_bw() +
    guides(fill = F) + xlab('') + ylab('Number of genes') +
    geom_boxplot(fill = 'white', width = .1) + axisFontSizes +
    scale_fill_manual(values = race.col)

counts.by.gender = ggplot(data = ind.counts, aes(x = GENDER, y = COUNT)) +
    geom_violin(aes(fill = GENDER), scale = 'width', width = .5) + theme_bw() +
    guides(fill = F) + xlab('') + ylab('Number of genes') +
    geom_boxplot(fill = 'white', width = .1) + axisFontSizes

counts.by.age = ggplot(data = ind.counts, aes(x = AGE, y = COUNT)) +
    geom_point(colour = 'dodgerblue3') + geom_smooth(colour = 'black') + theme_bw() +
    xlab('Age (years)') + ylab('Number of genes') + axisFontSizes

counts.by.bmi = ggplot(data = ind.counts, aes(x = BMI, y = COUNT)) +
    geom_point(colour = 'dodgerblue3') + geom_smooth(colour = 'black') + theme_bw() +
    xlab('BMI') + ylab('Number of genes') + axisFontSizes

counts.by.isch = ggplot(data = ind.counts, aes(x = ISCH, y = COUNT)) +
    geom_point(colour = 'dodgerblue3') + geom_smooth(colour = 'black') + theme_bw() +
    xlab('Ischemic time (minutes)') + ylab('Number of genes') + axisFontSizes

## load figure of rare variant enrichments produced by figure2a.count.enrichments.R 
load(paste0(dir, '/data/figure2a.count.enrichments.RData'))

count.ratio.plot.outliers.thresh = count.ratio.plot.outliers.thresh + axisFontSizes +
    scale_alpha_discrete(range = c(1, 0.4), name = 'Excluded individuals',
                         labels = c('None',
                                    expression('Outlier for' >= '50 genes'),
                                    expression('Outlier for' >= '30 genes'))) +
    theme(legend.text.align = 0)

outlier.bias = ggdraw() +
    draw_plot(ind.counts.plot,
              0/3, 0.7, 1/3, 0.3) +
    draw_plot(counts.by.race,
              1/3, 0.7, 1/3, 0.3) +
    draw_plot(counts.by.gender,
              2/3, 0.7, 1/3, 0.3) +
    draw_plot(counts.by.bmi,
              0/3, 0.4, 1/3, 0.3) +
    draw_plot(counts.by.age,
              1/3, 0.4, 1/3, 0.3) +
    draw_plot(counts.by.isch,
              2/3, 0.4, 1/3, 0.3) +
    draw_plot(count.ratio.plot.outliers.thresh, 0, 0, 1, 0.38) +
    draw_plot_label(c("a", "b", "c", "d", "e", "f","g"),
                    c(0, 1/3, 2/3, 0, 1/3, 2/3, 0),
                    c(1, 1, 1, 0.7, 0.7, 0.7, 0.37), size = 11)

pdf(paste0(dir, '/paper_figures/suppfig.number.outliers.by.covariates.pdf'), height = 7.5, width = 9, bg = 'white')

outlier.bias

dev.off()

# Test for significant effects of these covariates on number of outliers per individual
# RACE (compare White and Black)
race.test = wilcox.test(ind.counts$COUNT[ind.counts$RACE == 'White'], ind.counts$COUNT[ind.counts$RACE == 'Black'])
race.test
# W = 12209, p-value = 0.9281

# GENDER
gender.test = wilcox.test(COUNT ~ GENDER, data = ind.counts)
gender.test
# W = 22906, p-value = 0.5614

# AGE 
age.test = cor.test(ind.counts$COUNT, ind.counts$AGE, method = 's')
age.test
#S = 12846000, p-value = 0.03337
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1013435 

# BMI
bmi.test = cor.test(ind.counts$COUNT, ind.counts$BMI, method = 's')
bmi.test
#S = 14151000, p-value = 0.8339
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.0100142 

# Ischemic time
isch.test = cor.test(ind.counts$COUNT, ind.counts$ISCH, method = 's')
isch.test
#S = 11789000, p-value = 0.0002165
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1752725 
