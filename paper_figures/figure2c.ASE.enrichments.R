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
# assumes the necessary data.table is called data
get.outliers.controls.wrapper = function(dfsingle, dfmulti, thresholds, column = 'Z') {
    df = get.outliers.controls(dfsingle, dfmulti, thresholds[1], column)
    for (t in thresholds[2:length(thresholds)]) {
        df = rbind(df, get.outliers.controls(dfsingle, dfmulti, t, column))
    }
    return(df)
}

get.outliers.controls = function(dfsingle, dfmulti, thresh, column) {
    # get set of non-outliers with data    
    if (thresh > 0) {
        keep.single = dfsingle[,column] >= thresh
        keep.multi = dfmulti[,column] >= thresh
    } else {
        keep.single = dfsingle[,column] <= thresh
        keep.multi = dfmulti[,column] <= thresh
    }
    outsingle = apply(dfsingle[keep.single, c('GENE','INDS','TISSUE_ABBR')], 1, paste, collapse = "_")
    outmulti = apply(dfmulti[keep.multi, c('GENE','INDS')], 1, paste, collapse = "_")
    
    # restrict to genes that have at least one outlier (single- or multi-tissue)
    genes = unique(data[data$tissue_outlier %in% outsingle | data$outlier %in% outmulti, GENE_ID])
    plot.data = data[GENE_ID %in% genes, list(ASE, tissue_outlier, outlier, GENE_ID, TISSUE_ID)]
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
    plot.data$keep = FALSE
    plot.data$GENE_TISSUE = paste(plot.data$GENE_ID, plot.data$TISSUE_ID, sep= "_")
    non.single = unique(plot.data$GENE_TISSUE[plot.data$type == "Non-outlier"])
    single = unique(plot.data$GENE_TISSUE[plot.data$type == "Single-tissue outlier"])
    acceptable.pairs = intersect(non.single, single)
    # keep single-tissue outliers with a marked non-outlier
    plot.data$keep[plot.data$type == "Single-tissue outlier" & plot.data$GENE_TISSUE %in% acceptable.pairs] = TRUE

    non.multi = unique(plot.data$GENE_ID[plot.data$type == "Non-outlier"])
    multi = unique(plot.data$GENE_ID[plot.data$type == "Multi-tissue outlier"])
    acceptable.genes = intersect(non.multi, multi)
    # keep multi-tissue outliers with a matched non-outlier
    plot.data$keep[plot.data$type == "Multi-tissue outlier" & plot.data$GENE_ID %in% acceptable.genes] = TRUE
    # keep non-outliers that serve as controls for either single-tissue or multi-tissue outliers
    plot.data$keep[plot.data$type == "Non-outlier" & (plot.data$GENE_TISSUE %in% acceptable.pairs |
                                                      plot.data$GENE_ID %in% acceptable.genes)] = TRUE
    # remove all instances that aren't marked for keeping
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
# remove low mappability sites
ase.data = ase.data[ase.data$LOW_MAPABILITY==0,]
# add labeling sites near rare variants (chr 23 is X, but doesn't matter because removed chrX from ASE data)
nearvar = read.table(paste0(dir, '/preprocessing/ASEnearbyRare.02042017.txt'), header=T, stringsAsFactors=F)
rare.names = paste(nearvar$SUBJECT_ID, nearvar$CHR, nearvar$POS, sep='_')
ase.data[, near_rare:=(paste(SUBJECT_ID,CHR,POS,sep= "_") %in% rare.names)]
ase.data.with.rare = ase.data
# remove sites near rare variants
ase.data = ase.data[ase.data$near_rare==0,]
# add columns with ASE
ase.data$ASE = abs(ase.data$REF_RATIO - 0.5)

## Process outliers
outlierdir = paste0(dir,'/data/')
medz = read.picked(paste0(outlierdir, 'outliers_medz_nothreshold_picked.txt'))
tissuefiles = list.files(paste0(dir, '/data/singlez'), '*_picked.txt', full.names = T)
singletissues = lapply(tissuefiles, read.picked)
singlez = singletissues[[1]]
for (i in 2:length(singletissues)) {
    singlez = rbind(singlez, singletissues[[i]])
}
tissues = read.table(paste0(outlierdir, 'gtex_tissue_colors.txt'), header = T, stringsAsFactors = F, sep = "\t")
singlez$TISSUE_ABBR = tissues[match(singlez$TISSUE, tissues$tissue_site_detail_id), 'tissue_site_detail_abbr']

## get outliers (and matched non-outliers), for a variety of Z thresholds
data = ase.data
plot.data.all = get.outliers.controls.wrapper(singlez, medz, thresholds = c(-5:-1,1:5))

## test for significance between outliers and non outliers
# do this separately for single- and multi-tissue outliers
for (out in c("Single-tissue outlier", "Multi-tissue outlier")) {
    for (t in c(1:5)) {
        cat("Comparing", out, "to Non-outlier for threshold", t, ":\n")
        print(wilcox.test(plot.data.all$ASE[abs(plot.data.all$thresh) == t & plot.data.all$type == out],
                          plot.data.all$ASE[abs(plot.data.all$thresh) == t & plot.data.all$type == "Non-outlier"]))
    }
}
# all are significant p < 2.2e-16

# some switched to factors for plotting
plot.data.all$type = factor(plot.data.all$type, levels = c('Non-outlier', 'Single-tissue outlier', 'Multi-tissue outlier'),
    labels = c('Non-outlier', 'Single-tissue\noutlier', 'Multi-tissue\noutlier'))
plot.data.all$thresh = factor(abs(plot.data.all$thresh))

ase.colors = c('darkgrey','mediumorchid4','dodgerblue3')
names(ase.colors) = c('Non-outlier', 'Single-tissue\noutlier', 'Multi-tissue\noutlier')

ase.plot = ggplot(plot.data.all, aes(x = type, y = ASE, fill = type, alpha = thresh)) +
    geom_boxplot(outlier.colour = NA) +
    scale_alpha_discrete(range = c(0.3, 1),
                         name = 'Z-score threshold', breaks = c('1','2','3','4','5')) +
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

## FOR SUPPLMENTARY PANEL: look specifically at rare variants with the expectation that
## the ASE should be concordant with the observed outlier effect.
## For instance, for an underexpression outlier, we expect the rare variant to be underexpressed.

## read in rare variants associated with outliers and matched non-outliers
rarevars = fread(paste0(dir, '/features/byGene/10kb_genebody/all_rare_variants_SNPs.txt'))
colnames(rarevars) = c('SUBJECT_ID','GENE_ID','CHR','POS')
rarevars$GENE_ID= sapply(rarevars$GENE_ID, function (x) strsplit(x,".", fixed=T)[[1]][1])
rarevars$CHR = sub('chr', '', rarevars$CHR, fixed = T)
rarekeys = apply(rarevars, 1, paste, collapse = '_')
rarekeys = gsub(" ", "", rarekeys, fixed = T) # to deal with random spaces added in
ase.data.with.rare = ase.data.with.rare[paste(SUBJECT_ID,GENE_ID,CHR,POS,sep= "_") %in% rarekeys]

## add ase for rare variant ase matrix based on which allele is the minor allele
# let ase be the mirno allele ratio
alleles = fread(paste0(dir, '/features/variantBeds/AF_SNPs.bed'))
alleles = alleles[V4 <= 0.01 & V4 > 0, c(1,3,5), with = F]
setnames(alleles, c('CHR','POS','MINOR'))
alleles$CHR = as.character(alleles$CHR)
setkey(alleles, CHR, POS)
setkey(ase.data.with.rare, CHR, POS)
ase.data.with.rare = alleles[ase.data.with.rare]
ase.data.with.rare$ASE = ifelse(ase.data.with.rare$MINOR == 2, yes = 1 - ase.data.with.rare$REF_RATIO, no = ase.data.with.rare$REF_RATIO)

data = ase.data.with.rare
plot.data.rare = get.outliers.controls.wrapper(singlez, medz, thresholds = c(-3:-1,1:3))

## test for significance between outliers and non outliers
# do this separately for single- and multi-tissue outliers
for (out in c("Single-tissue outlier", "Multi-tissue outlier")) {
    for (t in c(-3:-1, 1:3)) {
        cat("Comparing", out, "to Non-outlier for threshold", t, ":\n")
        print(wilcox.test(plot.data.rare$ASE[plot.data.rare$thresh == t & plot.data.rare$type == out],
                          plot.data.rare$ASE[plot.data.rare$thresh == t & plot.data.rare$type == "Non-outlier"]))
    }
}

# some switched to factors for plotting
plot.data.rare$type = factor(plot.data.rare$type, levels = c('Non-outlier', 'Single-tissue outlier', 'Multi-tissue outlier'),
                             labels = c('Non-outlier', 'Single-tissue\noutlier', 'Multi-tissue\noutlier'))
plot.data.rare$thresh.abs = factor(abs(plot.data.rare$thresh))
plot.data.rare$class = ifelse(plot.data.rare$thresh > 0, yes = 'Overexpression', no = 'Underexpression')
plot.data.rare$class = factor(plot.data.rare$class, levels = c('Underexpression', 'Overexpression'))
plot.data.rare$thresh = factor(plot.data.rare$thresh)

ase.plot.rare = ggplot(plot.data.rare,
                        aes(x = type, y = ASE, fill = type, alpha = thresh.abs)) +
    facet_grid(~class) +
    geom_boxplot(outlier.colour = NA) +
    geom_hline(yintercept = 0.5, size = 0.75, linetype = 'dashed') +
    scale_alpha_discrete(range = c(0.3,1),
                         name = 'Z-score threshold', breaks = c(1:3), labels = c(-1:-3)) +
    theme_bw() + guides(alpha = guide_legend(override.aes = list(fill = 'darkgrey'),
                                             title.position = "top")) + 
    scale_fill_manual(values = ase.colors, guide = F) +
    xlab('') + ylab('Minor allele ratio') + 
    theme(axis.title = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 11),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.position = c(0.15, 0.92),
          legend.direction = 'horizontal',
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 11, face = 'bold'))
 
## Save main and supplemental plots
save(ase.plot, ase.plot.rare, ase.colors, file = paste0(dir, '/data/figure2c.ASE.enrichments.RData'))
