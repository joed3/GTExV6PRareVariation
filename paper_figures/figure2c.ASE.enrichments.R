#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# comparing ASE between outliers and non outliers
# use the best outlier per gene only

require(data.table)
require(ggplot2)

############################################
# FUNCTIONS

# reads in file and strips everything after the period in gene names
read.picked = function(filename) {
    df = read.table(filename, header=T, stringsAsFactors=F)
    df$GENE = sapply(df$GENE, function (x) strsplit(x,".", fixed=T)[[1]][1])
    return(df)
}

# returns data frame with data required for plotting single- and multi-tissue outliers and controls
# controls are genes that are either single or multi-tissue outliers for individuals that are not outliers in either class
# uses absolute values [to deal with Z-scores] of the given column
# ASSUMES ase.data exists
get.outliers.controls = function(dfsingle, dfmulti, thresh, column='Z') {
    # get set of non-outliers with data    
    keep.single = abs(dfsingle[,column]) >= thresh
    outsingle = apply(dfsingle[keep.single, c('GENE','INDS','TISSUE_ABBR')], 1, paste, collapse = "_")
    keep.multi = abs(dfmulti[,column]) >= thresh
    outmulti = apply(dfmulti[keep.multi, c('GENE','INDS')], 1, paste, collapse = "_")
    
    # restrict to genes that have at least one outlier (single- or multi-tissue)
    genes = unique(ase.data[ase.data$tissue_outlier %in% outsingle | ase.data$outlier %in% outmulti, GENE_ID])
    plot.data = ase.data[GENE_ID %in% genes, list(ASE, tissue_outlier, outlier, GENE_ID, TISSUE_ID)]
    plot.data[, single_outlier := plot.data$tissue_outlier %in% outsingle]
    plot.data[, multi_outlier := plot.data$outlier %in% outmulti]
    plot.data = as.data.frame(plot.data)
    plot.data$type = 'Non-outlier'
    plot.data$type[plot.data$single_outlier] = 'Single-tissue outlier'
    plot.data$type[plot.data$multi_outlier] = 'Multi-tissue outlier'
    # deal with ones that are both single and multi
    extra = plot.data[plot.data$single_outlier & plot.data$multi_outlier, ]
    extra$type = 'Single-tissue outlier'
    plot.data = rbind(plot.data, extra)
    plot.data$thresh = thresh

    # remove outliers and non-outliers that don't have a matched pair
    plot.data$keep = TRUE
    plot.data.non = plot.data[plot.data$type == "Non-outlier", ]
    plot.data.non$GENE_TISSUE = paste(plot.data.non$GENE_ID, plot.data.non$TISSUE_ID, sep = "_")
    plot.data.single = plot.data[plot.data$type == "Single-tissue outlier", ]
    plot.data.single$GENE_TISSUE = paste(plot.data.single$GENE_ID, plot.data.single$TISSUE_ID, sep= "_")
    plot.data.multi = plot.data[plot.data$type == "Multi-tissue outlier", ]
    # first, mark single-tissue outliers without a non-outlier (same gene, tissue)
    single.rm = plot.data.single$tissue_outlier[!(plot.data.single$GENE_TISSUE %in% plot.data.non$GENE_TISSUE)]
    plot.data$keep[plot.data$type == "Single-tissue outlier" & plot.data$tissue_outlier %in% single.rm] = FALSE
    # second, mark multi-tissue outliers without a non-outlier (same gene)
    multi.rm = plot.data.multi$outlier[!(plot.data.multi$GENE_ID %in% plot.data.non$GENE_ID)]
    plot.data$keep[plot.data$type == "Multi-tissue outlier" & plot.data$outlier %in% multi.rm] = FALSE
    # finally, mark non-outliers without a single-tissue outlier (same gene, tissue) or a multi-tissue outlier (same gene)
    non.rm = plot.data.non$tissue_outlier[!(plot.data.non$GENE_TISSUE %in% plot.data.single$GENE_TISSUE |
                                            plot.data.non$GENE_ID %in% plot.data.multi$GENE_ID)]
    plot.data$keep[plot.data$type == "Non-outlier" & plot.data$tissue_outlier %in% non.rm] = FALSE
    # remove all instances that got marked
    plot.data = plot.data[plot.data$keep, ]
    
    return(plot.data)
}

############################################

## Process ASE data with appropriate filters
# read in ASE data (and combine into a single data frame)
ase.data = fread(paste0('zcat ', dir, '/features/annotations/ASE/ase.v6.txt.gz'), sep='\t', header=T)
# remove chrX sites
ase.data = ase.data[ase.data$CHR %in% as.character(c(1:22)),]
# remove site for individuals with >=50 outliers (we exclude these both from the outlier set and the controls)
indiv.outlier.counts = read.table(paste0(dir, '/data/outliers_medz_picked_counts_per_ind.txt'),
    header = F, stringsAsFactors = F)
indivs.keep = indiv.outlier.counts[indiv.outlier.counts[,2] < 50, 1]
ase.data = ase.data[ase.data$SUBJECT_ID %in% indivs.keep,]
# only keep sites with at least 30 reads total
# and at least 5 ref and 5 nonref sites
ase.data = ase.data[ase.data$TOTAL_COUNT >= 30 & ase.data$REF_COUNT >=5 & ase.data$ALT_COUNT >=5,]
ase.data[, outlier := paste0(GENE_ID,'_',SUBJECT_ID)]
ase.data[, tissue_outlier := paste0(GENE_ID,'_',SUBJECT_ID,'_',TISSUE_ID)]
setkey(ase.data,outlier)
# remove biased sites
ase.data = ase.data[ase.data$MAPPING_BIAS_SIM==0,]
# add labeling sites near rare variants (chr 23 is X, but doesn't matter because removed chrX from ASE data)
nearvar = read.table(paste0(dir, '/preprocessing/ASEnearbyRare.09242015.txt'), header=T, stringsAsFactors=F)
rare.names = paste(paste0('GTEX-',nearvar$INDid), nearvar$chr_bps_siteAll_1, nearvar$chr_bps_siteAll_2, sep='_')
ase.data[, near_rare:=(paste(SUBJECT_ID,CHR,POS,sep= "_") %in% rare.names)]
# remove sites near rare variants
ase.data = ase.data[ase.data$near_rare==0,]
# remove low mappability sites
ase.data = ase.data[ase.data$LOW_MAPABILITY==0,]
# add columns with ASE
ase.data$ASE = abs(ase.data$REF_RATIO - 0.5)

## Process outliers
outlierdir = paste0(dir,'/data/')
medz = read.picked(paste0(outlierdir, 'outliers_medz_picked.txt'))
single = read.picked(paste0(outlierdir, 'outliers_singlez_picked.txt'))
tissues = read.table('gtex_tissue_colors.txt', header = T, stringsAsFactors = F, sep = "\t")
single$TISSUE_ABBR = tissues[match(single$TISSUE, tissues$tissue_site_detail_id), 'tissue_site_detail_abbr']

# get outliers (and matched non-outliers), for a variety of Z thresholds
plot.data.all = get.outliers.controls(single, medz, thresh = 5)
plot.data.all = rbind(plot.data.all, get.outliers.controls(single, medz, thresh = 4))
plot.data.all = rbind(plot.data.all, get.outliers.controls(single, medz, thresh = 3))
plot.data.all = rbind(plot.data.all, get.outliers.controls(single, medz, thresh = 2))

# test for significance between outliers and non outliers
# do this separately for single- and multi-tissue outliers
for (out in c("Single-tissue outlier", "Multi-tissue outlier")) {
    for (t in c(2:5)) {
        cat("Comparing", out, "to Non-outlier for threshold", t, ":\n")
        print(wilcox.test(plot.data.all$ASE[plot.data.all$thresh == t & plot.data.all$type == out],
                          plot.data.all$ASE[plot.data.all$thresh == t & plot.data.all$type == "Non-outlier"]))
    }
}

# some switched to factors for plotting
plot.data.all$type = factor(plot.data.all$type, levels = c('Non-outlier', 'Single-tissue outlier', 'Multi-tissue outlier'),
    labels = c('Non-outlier', 'Single-tissue\noutlier', 'Multi-tissue\noutlier'))
plot.data.all$thresh = factor(plot.data.all$thresh)

ase.colors = c('darkgrey','mediumorchid4','dodgerblue3')
names(ase.colors) = c('Non-outlier', 'Single-tissue\noutlier', 'Multi-tissue\noutlier')

ase.plot = ggplot(plot.data.all, aes(x = type, y = ASE, fill = type, alpha = thresh)) +
    geom_boxplot(outlier.colour = NA) +
    scale_alpha_discrete(range = c(0.3, 1),
                         name = 'Z-score threshold', breaks = c('2','3','4','5')) +
    theme_bw() + guides(alpha = guide_legend(override.aes = list(fill = 'darkgrey'),
                                             title.position = "top")) + 
    scale_fill_manual(values = ase.colors, guide = F) +
    xlab('') + ylab('ASE') + 
    theme(axis.title = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 11),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.position = c(0.3, 0.8),
          legend.direction = 'horizontal')

# Save workspace image
save(ase.plot, ase.colors, file = paste0(dir, '/data/figure2c.ASE.enrichments.RData'))

