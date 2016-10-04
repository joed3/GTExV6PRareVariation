#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

## Load required packages
require(ggplot2)
require(reshape2)
require(cowplot)

#------------- FUNCTIONS
trim.leading <- function (x)  sub("^\\s+", "", x)

#------------- MAIN

# Read in list of individuals with WGS
wgs = read.table('../preprocessing/gtex_2015-01-12_wgs_ids.148.txt', sep = '\t', header = F, stringsAsFactors = F)[, 1]
wgs.euro = scan('../preprocessing/gtex_2015-01-12_AF_ids.txt', what = character())

# Read in processed Sample ID and tissue file
meta = read.table(paste0(dir, '/preprocessing/gtex_2015-01-12_tissue_by_ind.txt'), header = T, sep = '\t', stringsAsFactors = F)
meta$INDid = trim.leading(meta$Id) 

# Get tissue names
tissues = unique(meta$Tissue)
tissues.clean = gsub('_', ' ', tissues)

# Get sample names
samples = unique(meta$INDid)
# Order samples by WGS or not (and, within wgs, euro or not)
samples.wgs = samples %in% wgs + 0
samples.wgs = samples.wgs + samples %in% wgs.euro
samples = samples[order(-samples.wgs)]
# clean up names
samples.clean = gsub('GTEX-', '', samples) 
wgs = gsub('GTEX-', '', wgs)
wgs.euro = gsub('GTEX-', '', wgs.euro)
wgs.noteuro = wgs[(!wgs %in% wgs.euro)]

# Make matrix of samples by tissues with expression data
exp.design = matrix(0, ncol = length(samples), nrow = length(tissues))
rownames(exp.design) = tissues.clean
colnames(exp.design) = samples.clean
col.counter = 1
for(i in samples){
	row.counter = 1
	ind.seq = meta[meta$Id == i, ]
	for(j in tissues){
		exp.design[row.counter, col.counter] = (j %in% ind.seq$Tissue) 
		row.counter = row.counter + 1
	}
	col.counter = col.counter + 1
}

# Cluster the matrix by rows
tissue.distances = as.dist(1 - cor(t(exp.design), method = 's'))
clustered.rows = hclust(tissue.distances)
exp.design = exp.design[clustered.rows$order, ]

# Cluster the matrix by columns
# First cluster the individuals with WGS
# Then cluster the individuals with only OMNI data
wgs.noteuro.distances = as.dist(1 - cor(exp.design[, colnames(exp.design) %in% wgs.noteuro], method = 's'))
clustered.wgs.noteuro = hclust(wgs.noteuro.distances)

wgs.euro.distances = as.dist(1 - cor(exp.design[, colnames(exp.design) %in% wgs.euro], method = 's'))
clustered.wgs.euro = hclust(wgs.euro.distances)

omni.distances = as.dist(1 - cor(exp.design[, !(colnames(exp.design) %in% c(wgs.euro, wgs.noteuro))], method = 's'))
clustered.omni = hclust(omni.distances)

exp.design = exp.design[, c(clustered.wgs.euro$order,
    length(clustered.wgs.euro$order) + clustered.wgs.noteuro$order,
    length(c(clustered.wgs.euro$order, clustered.wgs.noteuro$order)) + clustered.omni$order)]

# Melt the matrix for easy plotting with geom_tile
exp.design.melted = melt(exp.design)
colnames(exp.design.melted) = c('Tissue', 'Individual', 'Sampled')
exp.design.melted$Tissue = factor(as.character(exp.design.melted$Tissue), levels = rev(levels(exp.design.melted$Tissue)))
exp.design.melted$WGS = ifelse(exp.design.melted$Individual %in% c(wgs.noteuro, wgs.euro), 'WGS', 'OMNI')
exp.design.melted$WGS[exp.design.melted$Individual %in% wgs.euro] = 'WGS.EURO'

# Set alpha levels for those with and without WGS
colors = c('#FFC7B2','orangered','orangered4')
names(colors) = c('OMNI', 'WGS', 'WGS.EURO')

gtex.design = ggplot(data = exp.design.melted, aes(y = Tissue, x = Individual)) +
    geom_tile(aes(alpha = factor(Sampled), fill = factor(WGS)), colour = NA, size = 0.01) +
    scale_alpha_manual(values = c(0,1)) + scale_fill_manual(values = colors) +
    guides(fill = F, alpha = F) + 
    theme(axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 11))

# Calculate mean and median number of tissues per individual
avg.ntissues = mean(colSums(exp.design))
cat("average nubmer of tissues", avg.ntissues, "\n")

med.ntissues = median(colSums(exp.design))
cat("median number of tissues", med.ntissues, "\n")

# make figure showing MAFs of rare variants 
load(paste0(dir, '/data/euro.subpop.RData'))
# data generated in ../feature_construction/assess.euro.subpop.R

# combine SNVs and indels
euro.plotdata = rbind(indel.plotdata, snv.plotdata)

euro.maf = ggplot(euro.plotdata, aes(x = 100 * MAF, fill = Subpopulation)) + geom_histogram(binwidth = 0.1) +
    facet_grid(.~Subpopulation) + theme_bw() + guides(fill = FALSE) +
    ylab('Number of variants') + xlab("Minor allele frequency (%)") +
    scale_y_continuous(labels = scales::comma) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(size = 11))

# combined supplementary figure
combined.design.maf = ggdraw() +
    draw_plot(gtex.design, 0,0.6,1,0.4) +
    draw_plot(euro.maf, 0,0,1,0.6) +
    draw_plot_label(c('a','b'), c(0,0), c(1, 0.6), size = 11)


pdf(paste0(dir, '/paper_figures/suppfig.gtex.design.euro.maf.pdf'), height = 3.25, width = 7.2)

combined.design.maf

dev.off()
