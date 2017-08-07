#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# putting together the components of figure 2

library(ggplot2)
library(cowplot)

# load premade figure panels
load(paste0(dir, '/data/figure2a.count.enrichments.RData'))
load(paste0(dir, '/data/figure2b.threshold.enrichments.RData'))
load(paste0(dir, '/data/figure2c.ASE.enrichments.RData'))

fontsize = 11
fontsizes = theme(axis.text = element_text(size = fontsize),
                  axis.title = element_text(size = fontsize),
                  strip.text = element_text(size = fontsize + 1),
                  legend.text = element_text(size = fontsize - 1),
                  legend.title = element_text(size = fontsize - 1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

# make uniform format
count.ratio.plot.snv.indel = count.ratio.plot.snv.indel + fontsizes +
    scale_y_continuous(breaks = c(1, 1.5, 2)) +
    theme(axis.title.x = element_text(hjust = 0.95))
count.ratio.plot.sv = count.ratio.plot.sv + fontsizes +
    scale_y_continuous(breaks = c(1, 25, 50, 75)) +
    theme(axis.title = element_text(colour = 'white'))

props.plot.medz.singlez = props.plot.medz.singlez + fontsizes + theme(legend.position = c(0.75, 0.85))

ase.plot = ase.plot + fontsizes +
    theme(plot.margin = unit(c(5.5,5.5,-5,5.5), "pt"))

figure2.combined = ggdraw() +
    draw_plot(count.ratio.plot.sv, 0.62, 0.55, 0.3, 0.45) +
    draw_plot(count.ratio.plot.snv.indel, 0.08, 0.55, 0.55, 0.45) +
    draw_plot(props.plot.medz.singlez, 0, 0, 0.53, 0.54) +
    draw_plot(ase.plot, 0.53, 0, 0.47, 0.54) +
    draw_plot_label(c('a','b','c'), c(0.08, 0, 0.53), c(1, 0.55, 0.55), size = 15)

pdf(paste0(dir, '/paper_figures/figure2.pdf'), width = 8.5, height = 6)

figure2.combined

dev.off()
