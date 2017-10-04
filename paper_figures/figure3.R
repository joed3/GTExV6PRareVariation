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
effects.plot = effects.plot + axisFontSizes + theme(plot.margin = margin(5.5,5.5,-5.5,0),
                                                    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1))

x1 = 0.049
x2 = 0.041
figure3.combined = ggdraw() + 
    draw_plot(effects.plot, 0.27, 0, 0.73, .6) +
    draw_plot(outlier.cats.forest.scaled, 0.005, 0.5, 0.27, 0.5) +
    draw_plot(add.effects.plot.scaled, 0.005, 0, 0.27, 0.5) +
    draw_line(c(0.018, 0.018), c(0.57, 0.81)) + # SNV/indel
    draw_line(c(x1, x1), c(0.81, 0.99)) + # SV
    draw_text(c('SNV / indel', 'SV'), c(0.01, x2), c(0.69, 0.9), size = 11, angle = 90) +
    draw_line(c(x1, x1), c(0.07, 0.26)) + # conservation
    draw_line(c(x1, x1), c(0.3, 0.44)) + # regulation
    draw_text(c('Conservation', 'Regulation'), c(x2, x2), c(0.17, 0.37), size = 11, angle = 90) +
    draw_plot_label(c('a', 'b', 'c', 'd', 'e'), c(0, 0, 0.28, 0.28, 0.28), c(1, 0.5, 1, .61, .37), size = 15)

pdf(paste0(dir, '/paper_figures/figure3.raw.pdf'), height = 6.2, width = 12, bg = 'white')

figure3.combined

dev.off()
