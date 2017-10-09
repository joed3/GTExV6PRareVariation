#!/usr/bin/env Rscript

rm(list = ls())

dir = Sys.getenv('RAREVARDIR')
annodir = paste0(dir, '/features/annotations')
exacdir = Sys.getenv('EXAC_DIR')

# Load required packages
require(ggplot2)
require(reshape2)
require(gtable)
require(gridExtra)
require(RColorBrewer)
require(plyr)
require(dplyr)

#--------------- FUNCTIONS

# Function to generate odds ratios comparing one gene set to a range of other gene sets
gene.ors.calc = function(data, type, column.name, gene.lists, gene.lists.names){
	ors.out = matrix(, ncol = 4, nrow = length(gene.lists.names))
	table.out = c()
	expected.out = c()
	for(i in 1:length(gene.lists.names)){
		data[, gene.lists.names[i]] = ifelse(data$gene %in% gene.lists[[i]] | data$hgnc %in% gene.lists[[i]], 1, 0)
		test.table = table(data[, column.name], data[, gene.lists.names[i]])
		test.expected = data.frame(List = gene.lists.names[i], 
			Expected = sum(test.table[2, ]) * sum(test.table[, 2]) / sum(test.table), Observed = test.table[2, 2])
		expected.out = rbind(expected.out, test.expected)
		list.test = fisher.test(test.table)
		rownames(test.table) = NULL
		table.out = rbind(table.out, test.table)
		ors.out[i, ] = c(list.test$estimate, list.test$conf.int, list.test$p.value)
	}
	colnames(ors.out) = c('OR', 'LOW', 'HIGH', 'PVAL')
	ors.out = as.data.frame(ors.out)
	ors.out$TYPE = type
	ors.out$SET = gene.lists.names
	ors.out = ors.out[order(ors.out$OR), ]
	ors.out$SET = factor(ors.out$SET, levels = rev(ors.out$SET)) 
	table.out = as.data.frame(table.out)
	colnames(table.out) = c('Background', 'Disease gene')
	table.out$List = rep(gene.lists.names, each = 2)
	table.out[, 'Outlier status'] = rep(c('Background', 'Outlier'), length(gene.lists.names))
	table.out = table.out[, c('List', 'Outlier status', 'Background', 'Disease gene')]
	return(list(ors = ors.out, table = table.out, expected = expected.out))
}

# function to write tables to include in the supplement
write.gene.table = function(fname, genelist) {
    columns.order = c('gene', 'hgnc', 'indiv', 'ntissue', 'medz')
    column.names = c('Gene', 'HGNC', 'GTExID', 'NumberOfTissues', 'MedianZscore')
    
    list.out = medz[medz$outlier == 1 & medz$hgnc %in% genelist, columns.order]
    colnames(list.out) = column.names
    write.table(list.out, paste0(dir, '/data/', fname), sep = '\t', col.names = T, row.names = F, quote = 1:3)
}

#--------------- MAIN

#-- Global constants

# Define gene set names
gene.lists.names = c('eGene', 
	'GWAS', 
	'OMIM', 
	'Orphanet', 
	'ClinVar', 
	'ACMG',
	'Cardio', 
	'Cancer', 
	'DDG2P', 
	'DDG2P Both DD and IF', 
	'DDG2P Confirmed', 
	'DDG2P Possible', 
	'DDG2P Probable', 
	'DD Cancer')

# Gene types to keep: lincRNA and protein coding
types.to.keep = c('lincRNA', 'protein_coding')

# Outlier threshold
medz.threshold = 2

#-- Analysis

# Read in table of ENSG and HGNC ids from GTEx V6P GTF
ensg.hgnc = read.table(paste0(dir, '/data/gtex.v6p.hgnc.ensg.mapper.txt'), 
	sep = '\t', header = F, stringsAsFactors = F)
names(ensg.hgnc) = c('gene', 'hgnc')

# Read in list of genes with type info
# Restrict to protein coding and lincRNA genes
# Combine with HGNC symbols
ensg.types = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), 
	sep = '\t', header = F, stringsAsFactors = F)
colnames(ensg.types) = c('gene', 'type')
ensg.types = ensg.types %>% filter(type %in% types.to.keep)

gene.info = merge(ensg.types, ensg.hgnc) %>%
	mutate(pretty = gsub('\\.[0-9]+', '', gene))

# Read in list of outlier tested genes
medz = read.table(paste0(dir, '/data/outliers_medz_nothreshold_picked.txt'), 
	sep = '\t', header = T, stringsAsFactors = F)
names(medz) = c('gene', 'indiv', 'ntissue', 'medz')

# Add HGNC symbols to medz data frame
# Define outlier genes
medz = medz %>% merge(., gene.info, all.x = T) %>%
	mutate(outlier = (abs(medz) >= medz.threshold) + 0)

# Read in Metasoft summary data
# Filter for autosomal lincRNA and protein coding genes
# Add HGNC symbols
meta = read.table(paste0(dir, '/data/GTExReleaseV6PMetatissue.summary.txt'), 
	sep = '\t', header = T, stringsAsFactors = F)
names(meta) = tolower(names(meta))
meta = meta %>% filter(gene %in% gene.info$gene) %>%
	merge(., gene.info, all.x = T)

# Define top shared eQTL genes 
Ntop = sum(medz$outlier)
meta = meta[order(meta$pvalue), ]
meta$egene = 0
meta$egene[1:Ntop] = 1
meta.genes = meta$gene[meta$egene == 1]

# Read in GWAS dataset
gwas = read.delim(paste0(annodir, '/GWAS/gwas_catalog_v1.0-downloaded_2015-11-30.tsv'), 
	header = T, sep = '\t', stringsAsFactors = F)

# Get list of reported genes
reported.genes = unique(sort(unlist(strsplit(gwas$REPORTED.GENE.S., ' - |, '))))

# Read in OMIM dataset
omim.genes = read.table(paste0(annodir, '/OMIM/omim.genes.txt'), 
	sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Read in OrphaNet dataset
orpha.genes.ensg = read.table(paste0(annodir, '/OrphaNet/orphanet.genes.txt'), 
	sep = '\t', header = F, stringsAsFactors = F)[, 1]
orpha.genes = gene.info$hgnc[gene.info$pretty %in% orpha.genes.ensg]

# Read in Cardio and Cancer genes 
cardio.gold.genes = read.table(paste0(annodir, '/GeneListsOther/cardio.genes.gold.standard.csv'), 
	sep = ',', header = T, stringsAsFactors = F)[, 1]
cancer.gold.genes = read.table(paste0(annodir, '/GeneListsOther/cancer.genes.gold.standard.csv'), 
	sep = ',', header = T, stringsAsFactors = F)[, 1]

# Read in ACMG data
acmg = read.table(paste0(annodir, '/ACMG/completed_acmg.csv'), 
	sep = ',', header = T, stringsAsFactors = F)
acmg.genes = acmg$Gene

# Read in the ClinVar genes
clin = read.csv(paste0(annodir, '/GeneListsOther/gene_condition_source_id'), 
	sep = '\t', stringsAsFactors = F)
clin.genes = clin$GeneSymbol

# Read in the DDG2P gene set
ddg2p = read.table(paste0(annodir, '/DDG2P/DDG2P_2_8_2017.csv.gz'), sep = ',', header = T, stringsAsFactors = F)
ddg2p.categories = names(table(ddg2p$DDD.category))
ddg2p.genes = unique(ddg2p$gene.symbol)

# Make list of gene lists to test
gene.lists = list(meta.genes,
	reported.genes, 
	omim.genes, 
	orpha.genes,
	clin.genes,  
	acmg.genes, 
	cardio.gold.genes, 
	cancer.gold.genes, 
	ddg2p.genes)

# Calculate ORs for outlier and shared eQTL genes
medz.ors = gene.ors.calc(medz, 'Outlier', 'outlier', gene.lists, gene.lists.names)
meta.ors = gene.ors.calc(meta, 'eGene', 'egene', gene.lists[-1], gene.lists.names[-1])

# Write out table of counts
write.table(medz.ors$table, paste0(dir, '/data/medz.disease.genes.table.txt'), sep = '\t', row.names = F, col.names = T, quote = F)
write.table(meta.ors$table, paste0(dir, '/data/metasoft.disease.genes.table.txt'), sep = '\t', row.names = F, col.names = T, quote = F)

# Combine the OR dataframes for easy plotting
full.ors = rbind(medz.ors$ors[-which(medz.ors$ors == 'eGene'), ], meta.ors$ors)
full.ors$TYPE = factor(full.ors$TYPE, levels = c('eGene', 'Outlier'))

# Plot the enrichments for outlier genes
# Remove GTEx eGenes and DDG2P genes
bad.lists = c(gene.lists.names[1])

medz.ors.plot = ggplot(data = filter(medz.ors$ors, !(SET %in% bad.lists)), aes(x = SET, y = OR)) + 
	geom_abline(intercept = 1, slope = 0, colour = 'darkgrey') + 
	geom_pointrange(aes(x = SET, ymin = LOW, ymax = HIGH), colour = 'dodgerblue3', size = 0.3) + 
	theme_bw() + coord_flip() + ylab('Odds ratio') + xlab('') + theme(plot.title = element_text(hjust = .5))
medz.ors.plot

# Plot enrichments for outlier genes and multi-tissue eGenes
full.ors.plot.alt = ggplot(data = filter(full.ors, !(SET %in% bad.lists)), aes(x = SET, y = OR, shape = TYPE)) +
    geom_abline(intercept = 1, slope = 0, colour = 'darkgrey') + 
    coord_flip() + geom_pointrange(aes(x = SET, ymin = LOW, ymax = HIGH), colour = 'dodgerblue3', size = .3, position = position_dodge(width = 0.7)) +
    theme_bw() + ylab('Odds ratio') + xlab('') + 
    guides(shape = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position = c(.8, .9),
          legend.key = element_blank(),
          legend.text = element_text(size = 10),
          legend.background = element_blank(), 
          plot.title = element_text(hjust = .5))
full.ors.plot.alt

# Write out information for each disease list
write.gene.table('gwas.reported.outliers.txt', reported.genes)
write.gene.table('omim.outliers.txt', omim.genes)
write.gene.table('orpha.outliers.txt', orpha.genes)
write.gene.table('clinvar.outliers.txt', clin.genes)
write.gene.table('acmg.outliers.txt', acmg.genes)
write.gene.table('cardio.outliers.txt', cardio.gold.genes)
write.gene.table('cancer.outliers.txt', cancer.gold.genes)
write.gene.table('ddg2p.all.outliers.txt', ddg2p.genes)

# Save workspace image
save.image(paste0(dir, '/data/figure4c.gene.list.enrichments.RData'))

