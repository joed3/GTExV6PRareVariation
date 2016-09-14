#!/usr/bin/env Rscript
# Load packages

require(ggplot2)
require(reshape2)
require(plyr)
require(gplots)
require(data.table)
require(foreach)
require(doMC)
require(matrixStats)
require(grid)
require(gridExtra)

# Register backend for parallelization
registerDoMC(cores = 10)

#------------ FUNCTIONS

# Generate heatmap or correlation matrix using ggplot
gg.heat <- function(cor.mat, title, low = 'white', high = 'dodgerblue4'){
	cor.dist = as.dist(1 - cor.mat)
	cluster = hclust(cor.dist, method = 'complete')
	cor.ordered = cor.mat[cluster$order, cluster$order]
	cor.melted = melt(cor.mat)
	cor.melted$Var1 = factor(as.character(cor.melted$Var1), levels = colnames(cor.ordered))
	cor.melted$Var2 = factor(as.character(cor.melted$Var2), levels = colnames(cor.ordered))
	corr.plot = ggplot(data = cor.melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile() +
		scale_x_discrete(breaks = NULL) + xlab('') + ylab('') + scale_fill_gradient(low = low, high = high, name = 'Spearman rho') +
		ggtitle(title) + scale_y_discrete(limit = rev(levels(cor.melted$Var1)))
	corr.plot
	return(corr.plot)
}

# Function to generate PCA summary plots for gene expression data
gg.pca <- function(expression, title){
	exp.pca = svd(scale(t(expression)))
	per.exp = 100 * exp.pca$d^2 / sum(exp.pca$d^2)
	scree.plot = ggplot(data = as.data.frame(per.exp), aes(x = 1:length(tissues), y = per.exp)) + geom_point(colour = 'dodgerblue3') +
		ylab('Percent Variance Explained') + xlab('PC') + theme_bw() + ggtitle(title)
	exp.pca.df = data.frame(PC1 = exp.pca$u[, 1], PC2 = exp.pca$u[, 2], TISSUE = tissues)
	biplot = ggplot(data = exp.pca.df, aes(x = PC1, y = PC2, label = TISSUE)) + geom_point(aes(colour = TISSUE)) + theme_bw() +
		ggtitle(title) + xlab(paste('PC1 (', round(per.exp[1], 3), '%)', sep = '')) +
		ylab(paste('PC2 (', round(per.exp[2], 3), '%)', sep = '')) + scale_colour_manual(name = 'Tissue', values = gtex.colors) + guides(colour = F) +
		geom_text(check_overlap = T, nudge_y = .01, size = 3)
	return(list(Scree = scree.plot, Biplot = biplot))
}

# Function to arrange multiple plots and share legend in ggplot2
grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}



#-------------- MAIN

# Read in list of individuals and tissues
samples = scan('preprocessing/gtex_2015-01-12_individuals_all_normalized_samples.txt', what = character())
tissues = scan('preprocessing/gtex_2015-01-12_tissues_all_normalized_samples.txt', what = character())

# Read in the flat file 
rpkm = fread('preprocessing/gtex_2015-01-12_rpkm.txt')

# Gene ID as key
keys = c('Tissue', 'Gene')
setkeyv(rpkm, keys)

# Read in list of all expressed genes 
expressed.genes = read.table('preprocessing/gtex.expressed.genes.txt', sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Filter RPKM flat file for only these expressed genes
rpkm = rpkm[Gene %in% expressed.genes]

# Calculate mean gene expression levels
# Takes a few minutes
rpkm_mean = ddply(rpkm, .(Gene, Tissue), function(x) mean(unlist(x[1, 3:ncol(x)]), na.rm = T), .parallel = T)

# Unmelt the datasets to yield matrices of genes by tissues
rpkm_mean_mat = dcast(rpkm_mean, Gene ~ Tissue)

rownames(rpkm_mean_mat) = rpkm_mean_mat$Gene

rpkm_mean_mat = rpkm_mean_mat[, -1]

rpkm_mean_mat = as.matrix(rpkm_mean_mat)

# Calculate mean/median/sd across tissues
tissue_means = rowMeans(rpkm_mean_mat)
tissue_sds = rowSds(rpkm_mean_mat)
tissue_meds = apply(rpkm_mean_mat, 1, median, na.rm = T)

tissue_stats = data.frame(GENE = rownames(rpkm_mean_mat), MEAN = tissue_means, MEDIAN = tissue_meds, SD = tissue_sds)

# Write out tissue stats and raw matrices
write.table(tissue_stats, 'data/genes.rpkm.summary.stats.txt', quote = F, sep = '\t', col.names = T, row.names = F)
write.table(rpkm_mean_mat, 'data/genes.rpkm.mean.mat.txt', quote = F, sep = '\t', col.names = T, row.names = T)

