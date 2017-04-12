#!/usr/bin/env Rscript

# Load required packages
require(ggplot2)
require(reshape2)
require(cowplot)
require(plyr)

### Master directory
dir = Sys.getenv('RAREVARDIR')

#------------- FUNCTIONS
make.sharing.plot = function(df) {
    p = ggplot(data = df, aes(x = variable, y = Row)) +
        geom_tile(aes(fill = value), colour = 'black') + theme_bw() +
	scale_fill_gradient(low = 'white', high = 'navyblue', limits = c(0, 0.33), 
                            name = 'Replication\nproportion', na.value = "grey65") + 
	theme(axis.ticks = element_blank(),
              axis.text = element_text(size = 0),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 0),
              plot.margin = unit(c(-3,-3,-3,-3), "mm"),
              panel.border = element_blank())
    return(p)
}

make.combined.plot = function(vert, horiz, sharing, legend) {
    combined = ggdraw() +
        draw_plot(vert, 0,0.347,0.05,0.654) +
        draw_plot(horiz, 0.005,0.95,0.55,0.05) +
        draw_plot(sharing, 0.05,0.4,0.6,0.55) +
        draw_plot(legend, 0.65,0,0.35,1)
    return(combined)
}

make.point.plot = function(tissuesdf, vertical = TRUE){
    if (vertical) {
        p = ggplot(tissuesdf, aes(x = 1, y = -order, label = tissue_site_detail))
    } else {
        p = ggplot(tissuesdf, aes(x = order, y = 1))
    }
    p = p + geom_point(aes(colour = tissue_site_detail), size = 2.8) +
        scale_colour_manual(values = colors) + guides(colour = FALSE) + 
        theme(axis.line = element_line(linetype = 0),
              axis.ticks = element_blank(),
              axis.text = element_text(size = 0),
              axis.title = element_text(size = 0))
    return(p)
}

make.scatter = function(vals1, vals2, xlabel, ylabel, colour = "mediumorchid4") {
    cor.res = cor.test(vals1, vals2)
    corval = signif(cor.res$estimate, 3)
    
    minval = min(c(vals1, vals2)) - 0.01
    maxval = max(c(vals1, vals2)) + 0.01
    
    scatter = ggplot(data.frame(xval = vals1, yval = vals2), aes(xval, yval)) +
        geom_point(colour = colour, size = 0.5) +
        geom_abline(slope = 1, intercept = 0, colour = "darkgrey") + theme_bw() +
        xlab(xlabel) + ylab(ylabel) + ylim(c(minval, maxval)) + xlim(c(minval, maxval)) +
        annotate("text", x = minval + 0.2*(maxval-minval), y = maxval - 0.02*(maxval-minval),
                 label = paste("Pearson's r =", corval), size = 3.3) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 11))
    return(scatter)
}

make.scatter.gradient = function(df) {
    cor.res = cor.test(df$value, df$orig.value)
    corval = signif(cor.res$estimate, 3)
    
    minval = min(c(df$value, df$orig.value)) - 0.01
    maxval = max(c(df$value, df$orig.value)) + 0.01
    
    ggplot(df, aes(x = value, y = orig.value, colour = Row.n)) +
        geom_abline(slope = 1, intercept = 0, colour = "darkgrey") + geom_point(size = 0.5) +
        scale_colour_gradient(low = "darkorchid4", high = "plum", name = "Full discovery\nsample size") +
        ylab('Replication in all individuals') +
        xlab('Replication in discovery individuals') + theme_bw() +
        ylim(c(minval, maxval)) + xlim(c(minval, maxval)) +
        annotate("text", x = minval + 0.2*(maxval-minval), y = maxval - 0.02*(maxval-minval),
                 label = paste("Pearson's r =", corval), size = 3.3) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = c(0.8, 0.3),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 11),
              legend.title = element_text(size = 10),
              legend.background = element_blank()) +
         guides(colour = guide_colourbar(barwidth = 0.7))
}

# returns a list of single-tissue replication plots (ggplots)
# either of length 3 (1.subsetted initial replication, 2.replication in discovery individuals, 3. scatter plot of comparison)
# or of length 5 (1.initial replication, 2.replication in discovery inds, 3.replication in distinct inds, 4.comparison 1&2, 5.comparison 2&3)
# the number of plots depends on the input matrix
plot.replications = function(repli) {
    # replace - by _ (for LCLs to match the sharing tissue name)
    repli$tissue1 = sub("-", "_", repli$tissue1, fixed = TRUE)
    repli$tissue2 = sub("-", "_", repli$tissue2, fixed = TRUE)
    pairs1in2 = paste(repli$tissue1, repli$tissue2, sep = "_")
    pairs2in1 = paste(repli$tissue2, repli$tissue1, sep = "_")
    pairs.keep = c(pairs1in2, pairs2in1)
    sharing.subset = sharing
    sharing.subset$pairs = paste(sharing.subset$Row, sharing.subset$variable, sep = "_")
    sharing.subset$value[!(sharing.subset$pairs %in% pairs.keep)] = NA

    # subset of intial plot using all data
    sharing.subset.plot = make.sharing.plot(sharing.subset)
    # discovery and replication in the same individuals
    # build dataframe with values to use
    sameind = data.frame(pair = c(pairs1in2, pairs2in1),
                         value = c(repli$tis1.in.tis2.shared, repli$tis2.in.tis1.shared),
                         stringsAsFactors = F)
    sharing.sameind = sharing.subset
    sharing.sameind$value[!is.na(sharing.sameind$value)] = sameind$value[match(sharing.sameind$pairs[!is.na(sharing.sameind$value)], sameind$pair)]
    sharing.sameind.plot = make.sharing.plot(sharing.sameind)
    # comparing original vs replication in same inds
    scatter.df = sharing.sameind[!is.na(sharing.sameind$value),]
    scatter.df$Row.n = tcounts$nt[match(scatter.df$Row, tcounts$Tissue)]
    scatter.df$orig.value = sharing.subset$value[!is.na(sharing.subset$value)]
    scatter.orig = make.scatter.gradient(scatter.df)
    
    if (!("tis1.in.tis2.only" %in% colnames(repli))) {
        return(list(sharing.subset.plot, sharing.sameind.plot, scatter.orig))
    } else {
        # discovery and replication in distinct sets of individuals
        diffind = data.frame(pair = c(pairs1in2, pairs2in1),
                             value = c(repli$tis1.in.tis2.only, repli$tis2.in.tis1.only),
                             stringsAsFactors = F)
        sharing.diffind = sharing.subset
        sharing.diffind$value[!is.na(sharing.diffind$value)] = diffind$value[match(sharing.diffind$pairs[!is.na(sharing.diffind$value)], diffind$pair)]
        sharing.diffind.plot = make.sharing.plot(sharing.diffind)
        # comparing replication in same vs diff individuals
        scatter.diff = make.scatter(sameind$value, diffind$value,
                                    'Replication in discovery individuals', 'Replication in distinct individuals')
        
        return(list(sharing.subset.plot, sharing.sameind.plot, sharing.diffind.plot, scatter.orig, scatter.diff))
    }
}


#------------- MAIN

# Read in tissue colors and names
tissues = read.table('gtex_tissue_colors.txt', header = T, stringsAsFactors = F, sep = "\t")

# Read in outlier sharing data
sharing.raw = read.table(paste0(dir, '/data/fromXin/Fig1_outlierSharing.txt'), header = T, stringsAsFactors = F, sep = '\t')
sharing.raw = sharing.raw[, -2]

# Melt the dataframe
sharing = melt(sharing.raw)
sharing$variable = factor(as.character(sharing$variable), levels = colnames(sharing.raw)[-1])
sharing$Row = factor(sharing$Row, levels = rev(levels(sharing$variable)))

sharing$value[is.na(sharing$value)] = NA # turn NaN into NA
sharing$value[sharing$Row == sharing$variable] = 1

# For prettier plotting (wider range of colors)
sharing.mod = sharing
sharing.mod$value[sharing.mod$value >= 0.3] = 0.3

# Process tissues for plotting (this is where you set the order of the tissues)
tissues$tissue_site_detail_id = sub("-", "_", tissues$tissue_site_detail_id, fixed = TRUE)
tissues$tissue_site_detail_id = factor(tissues$tissue_site_detail_id, levels = levels(sharing.mod$variable))
tissues$order = as.numeric(tissues$tissue_site_detail_id)
tissues = tissues[!is.na(tissues$order),]
colors = paste0("#",tissues$tissue_color_hex)
names(colors) = tissues$tissue_site_detail

colors.vertical = make.point.plot(tissues)

colors.horizontal = make.point.plot(tissues, vertical = FALSE)

sharing.plot = make.sharing.plot(sharing.mod)

legend = colors.vertical + geom_text(hjust = 0, x = 1.02, size = 3.5) +
    theme(plot.margin = margin(-20,0,-20,-220))

combined = make.combined.plot(colors.vertical, colors.horizontal, sharing.plot, legend)

# Make the sharing plot
pdf(paste0(dir, '/paper_figures/figure1b.outlier.sharing.raw.pdf'), height = 7.2, width = 8.5)

print(combined)

dev.off()

## SUPPLEMENTAL

# Make adjustments for the fact that some tissues have higher individual overlaps
# first remake the same heatmap, but removing tissue pairs that aren't tested in the controlled replication analysis

# Load data generated by call_outliers/single_tissue_replication.R
load(paste0(dir, '/data/single_tissue_replication.RData'))

# get number of samples per tissue - needed for scatter plots
tcounts = read.table(paste0(dir, '/preprocessing/gtex_2015-01-12_tissue_by_ind.txt'), header = T, stringsAsFactors=F)
tcounts$Tissue = sub('-', '_', tcounts$Tissue, fixed = T)
tcounts = ddply(tcounts, .(Tissue), summarize, nt = length(Id))

# 70 shared only
plots70 = plot.replications(replications.shared70)

# 70 both shared and disjoint
plots70.same.diff = plot.replications(replications70)

save.image(paste0(dir, '/data/figure1b.outlier.sharing.RData'))

