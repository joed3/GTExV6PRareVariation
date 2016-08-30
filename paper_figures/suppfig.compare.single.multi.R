#!/usr/bin/env Rscript

library(ggplot2)

### Master directory
dir = Sys.getenv('RAREVARDIR')

# make supplementary figure that looks at the number of single-tissue outliers that overlap with multi-tissue ones
# and compares that to the sample size in that tissue

fontsize = 10
fontsizes = theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(size = fontsize),
                  axis.title = element_text(size = fontsize))

# read in tissue colors
colors = read.table(paste0(dir, '/data/gtex_tissue_colors.txt'), header = T, stringsAsFactors = F, sep = "\t")

tis.colors = paste0("#", colors$tissue_color_hex)
names(tis.colors) = colors$tissue_site_detail_id

load(paste0(dir, '/data/suppfig.single.multi.overlap.RData'))

pdf(paste0(dir, '/paper_figures/suppfig.single.multi.overlap.pdf'), height = 4.5, width = 4.5)

ggplot(samples.tissue.overlap, aes(x = nsamples, y = prop, colour = TISSUE)) + geom_point() +
    theme_bw() + ylab('Proportion of single-tissue outliers that are multi-tissue outliers') + xlab('Number of samples') +
    guides(colour = FALSE) + scale_colour_manual(values = tis.colors) + fontsizes +
    annotate("text", x = 110, y = 0.055, label = paste("Pearson r =", signif(correlation$estimate,3)), size = 3.5) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

dev.off()
