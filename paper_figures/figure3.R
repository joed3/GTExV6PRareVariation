#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

##Load required packages
require(ggplot2)
require(cowplot)

#------------- MAIN

# Load forest and effect size plots
load(paste0(dir, '/data/figure3a.rare.variant.class.enrichments.RData'))
load(paste0(dir, '/data/figure3b.feature.enrichments.RData'))
load(paste0(dir, '/data/figure3de.outlier.effect.size.RData'))

# Make font sizes uniform
fsize = 11
axisFontSizes = theme(axis.text.x = element_text(size = fsize),
                      axis.text.y = element_text(size = fsize),
                      axis.title = element_text(size = fsize),
                      strip.text = element_text(size = fsize),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())

outlier.cats.forest.scaled = outlier.cats.forest.scaled + axisFontSizes
add.effects.plot.scaled = add.effects.plot.scaled + axisFontSizes
effects.plot = effects.plot + axisFontSizes + theme(plot.margin = margin(5.5,5.5,-5.5,0))

# Make figure with 5 panels
# Leave space to add panel C from Xin
mygrey = "grey45"
figure3.combined = ggdraw() + 
  draw_plot(outlier.cats.forest.scaled, 0.04, .79, .49, .21) +
  draw_plot(add.effects.plot.scaled, .59, .79, .41, .21) +
  draw_plot(effects.plot, 0, 0, 1, .34) +
  draw_line(c(0.05, 0.05), c(0.825, 0.92), colour = mygrey) + # SNV/indel
  draw_line(c(0.16, 0.16), c(0.925, 0.99), colour = mygrey) + # SV
  draw_text(c('SNV /\nindel', 'SV'), c(0.02, 0.14), c(0.87, 0.96), size = 11, colour = mygrey) +
  draw_line(c(0.63, 0.63), c(0.83, 0.9), colour = mygrey) + # conservation
  draw_line(c(0.64, 0.64), c(0.92, 0.975), colour = mygrey) + # regulation
  draw_text(c('Conservation', 'Regulation'), c(0.58, 0.6), c(0.87, 0.95), size = 11, colour = mygrey) +
  draw_plot_label(c('a', 'b', 'c', 'd', 'e'), c(0, 0.53, 0, 0, 0), c(1, 1, .79, .34, .2), size = 15)

pdf(paste0(dir, '/paper_figures/figure3.raw.pdf'), height = 13.5, width = 12, bg = 'white')

figure3.combined

dev.off()
