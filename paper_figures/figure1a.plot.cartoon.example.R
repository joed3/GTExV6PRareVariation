#!/usr/bin/env Rscript

# created the carcass of figure 1. the rest is done in inkscape.

library(ggplot2)
library(reshape2)

dir = Sys.getenv('RAREVARDIR')

expr = read.table('expr.subset.by.genes.txt', header = T, stringsAsFactors = F)
examples = read.table('possible.txt', header = F, stringsAsFactors = F)
colnames(examples) = c('gene','ind','ntissue','medz')

# limit to tissue for which we have a pretty cartoon
tissues = c("Artery_Aorta","Artery_Coronary","Esophagus_Gastroesophageal_Junction","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Stomach","Whole_Blood")

# for each gene, individual pair print the tissues if they are all contained in this set
get.tis = function(gene, ind) {
    indname = sub('-','.', ind, fixed = T)
    subset = expr[expr$Gene == gene, c('Tissue','Gene',indname)]
    subset = subset[!is.na(subset[,3]), ]
    stopifnot(nrow(subset) == 5)
    if (sum(subset$Tissue %in% tissues) == 4) {
        print(subset)
    }
}

dummy = mapply(get.tis, examples$gene, examples$ind)

expr.melted = melt(expr, id.vars = c('Tissue','Gene'))
expr.melted = expr.melted[!is.na(expr.melted$value), ]

# picked one that would have a unique cartoon for each tissue
indtissues = c("Adipose_Subcutaneous","Liver","Lung","Stomach","Whole_Blood")
plot.data = expr.melted[expr.melted$Gene == "ENSG00000198610.6" & expr.melted$Tissue %in% indtissues, ]

# for that gene, get the medz for all individual (with at least 5 tissues)
gene.expr = expr[expr$Gene == "ENSG00000198610.6", -c(1,2)]
gene.medz = apply(gene.expr, 2, function(x) ifelse(sum(!is.na(x)) >= 5, yes = median(x, na.rm = TRUE), no = NA))
gene.medz = gene.medz[!is.na(gene.medz)]
gene.medz.df = data.frame(Tissue = "Median", Gene = "ENSG00000198610.6", variable = names(gene.medz), value = gene.medz)

plot.data = rbind(plot.data, gene.medz.df)
plot.data$Tissue = factor(plot.data$Tissue, levels = c("Adipose_Subcutaneous","Liver","Stomach","Whole_Blood","Lung","Median"))

#cols = c("#FF6600","#AABB66","#99FF00","#FFDD99","#FF00BB")
#names(cols) = indtissues

pdf(paste0(dir, '/paper_figures/figure1a.cartoon.draft.pdf'), height = 7, width = 4.5)

ggplot(plot.data, aes(x = value)) +
    geom_histogram(binwidth = 0.15, colour = "white", fill = "darkgrey") + xlab('Z-score') + ylab('') +
    facet_grid(Tissue~., scales = "free") + theme_classic() + guides(fill = FALSE) + xlim(c(-5,5)) +
    geom_vline(xintercept = c(-0.9624364,-0.3787555,4.6037353,3.2303840,3.8570197), size = 1.1) +
    scale_fill_manual(values = cols) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 13))

dev.off()
