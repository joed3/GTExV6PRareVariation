#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

##Load required packages
require(ggplot2)
require(cowplot)

# Load ASE, gene enrichment and UK10K data
load(paste0(dir, '/data/figure4a.uk10k.RData'))
load(paste0(dir, '/data/figure4b.exac.enrichments.RData'))
load(paste0(dir, '/data/figure4c.gene.list.enrichments.RData'))

# Common axis font sizes
fsize = 11
axisFontSizes = theme(axis.text.y = element_text(size = fsize),
                      axis.text.x = element_text(size = fsize),
                      axis.title = element_text(size = fsize),
                      plot.title = element_text(size = fsize + 1),
                      legend.text = element_text(size = fsize - 1),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())

# Change X-axis labels for UK10K plot
uk10k.wleg = uk10k.wleg + scale_x_discrete(labels = c('0', '1', '2', '3', '4', '> 4')) +
    axisFontSizes + 
    theme(axis.text.x = element_text(angle = 0, hjust = .5),
          legend.position = c(0.5,0.75),
          legend.background = element_blank(),
          legend.key.size = unit(1, "lines"))

# Fix font sizes for alternative ExAC figure
ptv.odds.plot = ptv.odds.plot + axisFontSizes

# Use common font sizes for the disease gene enrichments
medz.ors.plot.alt = medz.ors.plot.alt + axisFontSizes

# Make figure with 3 panels: 1. UK10K, 2. ExAC, and 3. Disease gene enrichments
figure4.combined = ggdraw() + 
  draw_plot(uk10k.wleg, 0, 0, 0.32, 0.95) +
  draw_plot(ptv.odds.plot, 0.32, 0, 0.33, 0.95) +
  draw_plot(medz.ors.plot.alt, 0.65, 0, 0.35, 0.95) +
  draw_plot_label(c('a', 'b', 'c'), c(0, 0.33, 0.66), c(1, 1, 1), size = 15)
figure4.combined

pdf(paste0(dir, '/paper_figures/figure4.pdf'), width = 10, height = 3.375)

figure4.combined

dev.off()
