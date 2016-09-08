#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(ggplot2)

# Allele frequency colors
af.colors = c('darkgrey', 'dodgerblue3')
names(af.colors) = c('Non-outlier', 'Outlier')

# Load Xin's data for this figure
data = read.table(paste0(dir, '/data/fromXin/Fig4a_outlier_AF_uk10k.txt'), sep = '\t', header = T)

# modify factor labels
data$group = factor(data$group, levels = c('non-outlier','outlier'), labels = c('Non-outlier','Outlier'))

# Emulating the ExAC style
nind = 3781
#data$ACSimple = data$uk10k_ALTc
data$ACSimple = mapply(min, data$uk10k_ALTc, nind*2 - data$uk10k_ALTc)
data$ACSimple = ifelse(data$ACSimple > 4, 'AC>4', paste('AC=', data$ACSimple, sep = ''))
uk10k.table = table(data$ACSimple, data$group)
uk10k.table[, 1] = uk10k.table[, 1] / colSums(uk10k.table)[1]
uk10k.table[, 2] = uk10k.table[, 2] / colSums(uk10k.table)[2]
uk10k.df = as.data.frame(uk10k.table)
levels(uk10k.df$Var1) = c('0', '1', '2', '3', '4', '> 4')

uk10k.wleg = ggplot(data = uk10k.df, aes(x = Var1, y = Freq, fill = Var2)) +
	geom_bar(stat = 'identity', position = 'dodge', colour = 'black') + theme_bw() + xlab('Minor allele count in UK10K') +
	ylab('Proportion of promoter variants') + scale_fill_manual(values = af.colors) + theme(legend.title=element_blank()) +
	theme(legend.key = element_blank())
uk10k.wleg

# Use Wilcoxon rank sum test to obtain significance level
data$uk10k_MAF = mapply(min, data$uk10k_AF, 1-data$uk10k_AF)
uk10k.sig.test = wilcox.test(data$uk10k_MAF[data$group == 'Outlier'], data$uk10k_MAF[data$group == 'Non-outlier'])
uk10k.sig.test
#W = 30710, p-value = 0.001266

# Save workspace image
save.image(paste0(dir, '/data/figure4a.uk10k.RData'))
