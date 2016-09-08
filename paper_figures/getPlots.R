# ===================================================================================================================
# This is a script with necessary functions for some of the scripts for drawing figures
# ===================================================================================================================

# Recall required packages
library(ggplot2)
library(plotrix)

# 3D scatter plot
getPostProbs <- function(x, y, z) {
  plot(x, y, col = color.scale(1.5*z,c(0.75,0.1),c(0.75,.6),c(0.75,1)),
       xlab="Gene expression, |Median Z-score|", 
       ylab="Genomic annotations, P(FR | G)", 
       cex.axis=1.3, cex.lab=1.6, pch=16) # x & y
  
  rbPal = colorRampPalette(c('gray','dodgerblue')) # z
  Col = rbPal(20)[cut(1.5*z,breaks=20)]
  legend_image = as.raster(matrix(rbPal(20), nrow=1))
  rasterImage(legend_image,8,0.65,11,0.685)
  
  text(8.05,0.63,"0",cex=1.2)
  text(11.05,0.63,"1",cex=1.2)
  text(9.5,0.71,"P(FR | G, E)",cex=1.2)
  text(9.5,0.75,"RIVER",cex=1.2)
}

# scatter plot with examples
getScatterExmpls <- function(x, y, color) {
  for (i in 1:length(x)) {
    points(x[i], y[i], col = color, pch=19)
  }
}

# histograms of different groups based on relative frequencies
getHistExmpls <- function(data, bins) {
  ggplot(data, aes(x, fill=group, colour=group)) +
    geom_histogram(aes(y=2*(..density..)/sum(..density..)), 
                   breaks=bins, alpha=0.6, position="identity",lwd=0.2) +
    scale_color_manual(values=c("gray","darkorange2","firebrick"),
                       name="", breaks=c("A","B","C"),
                       labels=c("All","Pathogenic","Regulatory+Pathogenic"))+ 
    scale_fill_manual(values=c("gray","darkorange2","firebrick"),
                      name="", breaks=c("A","B","C"),
                      labels=c("All","Pathogenic","Regulatory+Pathogenic"))+ 
    theme_bw() + xlab("") + ylab("Relative frequency") + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size=8),
          axis.title.y = element_text(size=10),
          axis.text.x = element_blank(),
          legend.text = element_text(size=10))
}

# histograms of different groups based on relative frequencies (log10-transformed)
getHistExmpls2 <- function(data, xlabel, min, max) {
ggplot(dat, aes(x, fill=group, colour=group)) +
  geom_histogram(aes(y=2*(..density..)/sum(..density..)), 
                 breaks=seq(min, max, (max-min)/50), 
                 alpha=0.45,position="identity",lwd=0.2) +
  scale_color_manual(values=c("gray","darkorange2"),
                     name="",
                     breaks=c("A","B"),
                     labels=c("All","Pathogenic"))+
  scale_fill_manual(values=c("gray","darkorange2"),
                    name="",
                    breaks=c("A","B"),
                    labels=c("All","Pathogenic"))+
  theme_bw() + xlab(xlabel) + ylab("Relative frequency") +
  # annotate("text", x=-4, y=.15, label= "P < 1.7x10-3", size=4) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size=11))
}

# Multiple boxplots with groups having gradient colors 
getMultipleBoxPlotsGradients <- function(data) {
  ggplot(data, aes(x = type, y = value, fill = type, alpha = thrd)) +
    geom_hline(yintercept=0,color="gray") +
    geom_boxplot(outlier.color=NA)+
    scale_alpha_discrete(range = c(0.3, 1),
                         name = 'Z-score threshold', 
                         breaks = c('2','3','4')) +
    theme_bw() + 
    guides(alpha = guide_legend(title.theme = element_text(size=14,angle=0),
                                override.aes = list(fill = 'darkgrey'),title.position = "top")) + 
    scale_fill_manual(values = colors, guide = F) +
    xlab('') + ylab('Rank correlation') + 
    theme(axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.position = c(.2, .9),
          legend.text=element_text(size=12),
          legend.direction = 'horizontal')
}

# Forest plots with grouped colors
getForestPlots <- function(data, Annotation) {
  ggplot(data, aes(x=x,y=y), las=2)+
    geom_point(aes(colour=Annotation), cex=3)+
    geom_errorbarh(aes(xmin=xlow, xmax=xup, colour=Annotation), height=0.0, lwd=1)+
    theme_bw()+ xlab("Log odds ratio") + ylab("Genomic annotations")+
    scale_y_discrete(breaks=seq(1,nrow(data),1), labels=data$labels)+
    geom_vline(xintercept=0, color="gray", linetype="longdash") + # add a reference line
    scale_x_continuous(limits = c(-.3,.5))+
    scale_colour_manual(values=cbPalette)+
    theme(axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=9),
          axis.title.y = element_text(size=14),
          axis.text.y  = element_text(size=9),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top", 
          legend.title = element_blank(),
          legend.direction = "horizontal",
          legend.text = element_text(size=11),
          legend.key.size = unit(1.5, "lines"),
          legend.key = element_blank(),
          plot.margin = unit(c(1,2,1,1), "lines")) 
}

