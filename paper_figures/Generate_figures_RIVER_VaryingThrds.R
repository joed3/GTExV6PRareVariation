# ===================================================================================================================
# This is to generate figures for EDF9a
#
# input:
#       1. RIVER_VaryingThrds.out.RData
#          load all data used for running main_RIVER_VaryingThrds.R
#       2. RIVER_thrdx.out.RData (x: 1,1p5,2,2p5)
#          load all data used for running main_RIVER_VaryingThrds.R in each threshold case
#
# output:
#       1. auc_VaryingThrds_ROCcurve_x.pdf (x: 1,1p5,2,2p5)
#          AUC figures like figure 5b for each threshold case
#
# Note that the final figure was modified by combining four auc_VaryingThrds_ROCcurve_x.pdf figures via inkscape for visualization purpose
#
# ===================================================================================================================

# Master directory
dir = Sys.getenv('RAREVARDIR')

# Recall required functions
library(dplyr)
library(data.table)
library(ggplot2)
library(pROC)


load(file = paste("/data/RIVER_VaryingThrds.out.RData",sep=""))

## AUC (ROC curves)
thrds_table = data.frame(values=c(1,1.5,2,2.5),labels=c("1","1p5","2","2p5"))

# edge_type = c(2,1,4,3)
for (level in 1:nrow(thrds_table)) {
  pdf(paste("/data/auc_VaryingThrds_ROCcurve_",thrds_table[level,"labels"],".pdf",sep=""),useDingbats=FALSE, width=7, height=7)
  
  par(mar=c(6.1, 6.1, 4.1, 4.1))
  plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate", cex.axis=1.3, cex.lab=1.6)
  abline(0, 1, col="gray", lwd=2) 
  
  load(file = paste("/data/RIVER_thrd",as.character(thrds_table[level,"labels"]), ".out.RData",sep=""))
  lines(1-int.em.roc$specificities,int.em.roc$sensitivities, type="s", col='dodgerblue', lwd=4)
  lines(1-g.roc$specificities,g.roc$sensitivities, type="s", col='mediumpurple', lwd=4)
  roc.test(int.em.roc, g.roc)$p.value[1]
  # text(0.5, 0.95, bquote(Median ~ "|" ~ Z-score ~ "|" ~ ">=" ~ .(thrds_table[level,"values"])))
  text(0.3, 0.95, expression("Median |Z-score|" >= 1), cex=1.3)
  text(0.2, 0.88, expression(paste(italic("P"), " = ", "1.02 x ", 10^{-9}, sep="")), cex=1.3)
  legend(0.45, 0.2, c(paste("RIVER (AUC = ",round(int.em.roc$auc[1],3),")",sep=""),
                      paste("Genomic annotation\nmodel (AUC = ",round(g.roc$auc[1],3),")",sep="")), 
         lty=c(1,1), lwd=c(4,4), col=c("dodgerblue","mediumpurple"), cex=1.2, pt.cex=1.4, bty="n")
  dev.off()
}
