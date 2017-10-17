#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')
rpkmfile = 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct' # set to path of this downloaded file

## calculate median rpkm across tissue, including only tissues in which the gene is expressed
## really, it's the median across tissues of the median across individuals

library(data.table)
library(scales)

######################################################
## FUNCTIONS
clean.tissue.name = function(tissue) {
    tissue = gsub(' - ', '_', tissue, fixed = T)
    tissue = gsub(' ', '_', tissue, fixed = T)
    tissue = gsub('(', '', tissue, fixed = T)
    tissue = gsub(')', '', tissue, fixed = T)
    return(tissue)
}

######################################################

## first get the set of tissues each gene is expressed in
outfile = paste0(dir, '/preprocessing/gtex_2015-01-12_expressed_genes_by_tissue.txt')
command = paste0('cut -f1,2 ', dir, '/preprocessing/gtex_2015-01-12_normalized_expression.txt > ', outfile)
system(command)

## read in the median rpkms across individuals and the file created above
expressed = fread(outfile)
setkey(expressed, Gene)
rpkm = fread(rpkmfile, skip = 1)
## clean up column names of rpkm
tissuenames = colnames(rpkm)[3:ncol(rpkm)]
setnames(rpkm, tissuenames, clean.tissue.name(tissuenames))

## take medians across tissues, restricting to tissues with sufficient expression
expr.tissues = unique(expressed$Tissue)
stopifnot(sum(expr.tissues %in% colnames(rpkm)) == length(expr.tissues))

medians = data.frame(Gene = rpkm$Name, Median = NA, stringsAsFactors = F)
for (i in 1:nrow(medians)) {
    gene = rpkm[i, Name]
    stopifnot(medians[i, 'Gene'] == gene)
    tissues = unlist(expressed[gene, 'Tissue'])
    if (is.na(tissues[1])) {
        medians[i, 'Median'] = NA
    } else {
        medians[i, 'Median'] = median(unlist(rpkm[i, tissues, with = F]))
    }
}

write.table(medians, paste0(dir, '/preprocessing/gtex_2015-01-12_median_rpkm.txt'), col.names = T, row.names = F, quote = F, sep = '\t')
