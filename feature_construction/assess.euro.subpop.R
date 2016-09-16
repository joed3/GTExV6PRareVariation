#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# check the allele frequencies in european subpopulations of rare variants included in our analyses

library(data.table)
library(ggplot2)

# input files (use variant list created for the Battle lab: includes all variants within 10 kb of all genes)
snv = fread(paste0(dir, '/features/byGene/10kb/all_rare_variants_SNPs.txt'), header = F)
indel = fread(paste0(dir, '/features/byGene/10kb/all_rare_variants_indels.txt'), header = F)
setnames(snv, c('Ind','Gene','Chr','Pos'))
setnames(indel, c('Ind','Gene','Chr','Pos'))
setkey(snv, Gene)
setkey(indel, Gene)

# restrict to protein-coding and lincRNA genes
genetypes = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), stringsAsFactors = F, header = F)
genesToKeep = genetypes[genetypes[,2] %in% c('protein_coding','lincRNA'), 1]
snv = snv[Gene %in% genesToKeep, ]
indel = indel[Gene %in% genesToKeep, ]

# get unique set of variants of each type
snv = snv[, .(Chr, Pos)]
indel = indel[, .(Chr, Pos)]
setkey(snv) # reset key to all columns: allows use of unique
setkey(indel)
snv = unique(snv)
indel  = unique(indel)

# for each subpopulation, read in the appropriate AF file
popdir = paste0(dir, '/features/variantBeds/1KG/')
snv.af.list = list()
indel.af.list = list()
for (subpop in c('CEU','FIN','GBR','IBS','TSI')) {
    snv.af = fread(paste0(popdir, subpop, '.SNPs.AF.txt'), select = c(1,2,5))
    setnames(snv.af, c('CHROM','POS','AF'))
    indel.af = fread(paste0(popdir, subpop, '.indels.AF.txt'), select = c(1,2,5))
    setnames(indel.af, c('CHROM','POS','AF'))
    # add MAF
    snv.af = snv.af[, .(variant = paste0("chr", CHROM, "_", POS), MAF = min(AF, 1-AF)), by = .(CHROM, POS)]
    indel.af = indel.af[, .(variant = paste0("chr", CHROM, "_", POS), MAF = min(AF, 1-AF)), by = .(CHROM, POS)]
    # add vector of MAFs for variants actually in our analyses
    snv.af.subset = snv.af[variant %in% paste0(snv$Chr, "_", snv$Pos), MAF]
    snv.af.list[[subpop]] = c(snv.af.subset, rep(0, nrow(snv) - length(snv.af.subset)))
    indel.af.subset = indel.af[variant %in% paste0(indel$Chr, "_", indel$Pos), MAF]
    indel.af.list[[subpop]] = c(indel.af.subset, rep(0, nrow(indel) - length(indel.af.subset)))
}

# make data frame that ggplot can use
snv.plotdata = data.frame(Subpopulation = character(), MAF = numeric(), stringsAsFactors = F)
indel.plotdata = data.frame(Subpopulation = character(), MAF = numeric(), stringsAsFactors = F)
for (subpop in names(snv.af.list)) {
    snv.plotdata = rbind(snv.plotdata, data.frame(Subpopulation = subpop, MAF = snv.af.list[[subpop]], stringsAsFactors = F))
    indel.plotdata = rbind(indel.plotdata, data.frame(Subpopulation = subpop, MAF = indel.af.list[[subpop]], stringsAsFactors = F))
}

# actual plotting
pdf(paste0(dir, '/shared_figures/euro.maf.histogram.pdf'), height = 4, width = 10)

for (plotdata in list(snv.plotdata, indel.plotdata)) {
    print(ggplot(plotdata, aes(MAF, fill = Subpopulation)) + geom_histogram(binwidth = 0.001) + facet_grid(.~ Subpopulation) +
        theme_bw() + guides(fill = FALSE))
}

dev.off()

pdf(paste0(dir, '/shared_figures/euro.maf.density.pdf'), height = 5, width = 5)

for (plotdata in list(snv.plotdata, indel.plotdata)) {
    print(ggplot(plotdata, aes(MAF, fill = Subpopulation)) + geom_density(adjust = 10, alpha = 0.3) +
        theme_bw())
}

dev.off()

# compare MAFs between all pairs of subpopulations
subpops = c('CEU','FIN','GBR','IBS','TSI')
for (i in 1:4) {
    for (j in (i+1):5) {
        subpop1 = subpops[i]
        subpop2 = subpops[j]
        cat(subpop1,subpop2,'\n')
        cat('SNV\n')
        print(wilcox.test(snv.af.list[[subpop1]], snv.af.list[[subpop2]]))
        cat('indel\n')
        print(wilcox.test(indel.af.list[[subpop1]], indel.af.list[[subpop2]]))
    }
}

save(snv.plotdata, indel.plotdata, file = paste0(dir, '/data/euro.subpop.RData'))
