#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')
annodir = paste0(dir, '/features/annotations')
exacdir = Sys.getenv('EXAC_DIR')

# Load required packages
require(ggplot2)
require(RColorBrewer)

#--------------- FUNCTIONS

# Function to generate odds ratios comparing one gene set to a range of other gene sets
gene.ors.calc = function(data, type, column.name, gene.lists, gene.lists.names){
	ors.out = matrix(, ncol = 4, nrow = length(gene.lists.names))
	table.out = c()
	expected.out = c()
	for(i in 1:length(gene.lists.names)){
		data[, gene.lists.names[i]] = ifelse(data$GENE %in% gene.lists[[i]] | data$HGNC %in% gene.lists[[i]], 1, 0)
		test.table = table(data[, column.name], data[, gene.lists.names[i]])
		test.expected = data.frame(List = gene.lists.names[i], Expected = sum(test.table[2, ]) * sum(test.table[, 2]) / sum(test.table), Observed = test.table[2, 2])
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
    columns.order = c('GENE', 'HGNC', 'INDS', 'DFS', 'Z')
    column.names = c('Gene', 'HGNC', 'GTExID', 'NumberOfTissues', 'MedianZscore')
    
    list.out = medz[medz$OUTLIER == 1 & medz$HGNC %in% genelist, columns.order]
    colnames(list.out) = column.names
    write.table(list.out, paste0(dir, '/data/', fname), sep = '\t', col.names = T, row.names = F, quote = 1:3)
}

#--------------- MAIN

# Define gene set colors
gene.lists.names = c('eGene', 
	'GWAS', 
	'OMIM', 
	'Orphanet', 
	'ClinVar', 
	'ACMG',
	'Cardio', 
	'Cancer', 
	'LoF intolerant')

colors = brewer.pal(12, name = 'Paired')
colors = c(colors[c(1, 3:6, 12, 7:10)], 'orangered3')
colors = c(colors[1:5], 'red3', colors[6:length(colors)])

names(colors) = gene.lists.names

# Read in table of ENSG and HGNC ids from GTEx V6P GTF
ensg.hgnc = read.table(paste0(dir, '/data/gtex.v6p.hgnc.ensg.mapper.txt'), sep = '\t', header = F, stringsAsFactors = F, row.names = 1)

# Read in list of genes with type info
# Restrict to protein coding and lincRNA genes
ensg.types = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), sep = '\t', header = F, stringsAsFactors = F)
colnames(ensg.types) = c('GENE', 'TYPE')
ensg.types$HGNC = ensg.hgnc[ensg.types[, 1], 1]
rownames(ensg.types) = ensg.types$GENE
types.to.keep = c('lincRNA', 'protein_coding')
ensg.types = ensg.types[ensg.types$TYPE %in% types.to.keep, ]
ensg.types$PRETTY = gsub('\\.*', '', ensg.types$GENE)

# Read in list of outlier tested genes
medz = read.table(paste0(dir, '/data/outliers_medz_nothreshold_picked.txt'), sep = '\t', header = T, stringsAsFactors = F)

# Add HGNC symbols to medz data frame
medz$HGNC = ensg.hgnc[medz$GENE, 1]

# Define outlier genes
medz.threshold = 2
medz$OUTLIER = ifelse(abs(medz$Z) >= medz.threshold, 1, 0)

# Read in Metasoft summary data
meta = read.table(paste0(dir, '/data/GTExReleaseV6PMetasoft.summary.txt'), sep = '\t', header = T, stringsAsFactors = F)
colnames(meta)[1] = 'GENE'
meta = meta[meta$GENE %in% ensg.types$GENE, ]

# Add HGNC symbols to meta data frame
meta$HGNC = ensg.hgnc[meta$GENE, 1]

# Define top shared eQTL genes 
Ntop = sum(medz$OUTLIER)
meta = meta[order(meta$Pvalue), ]
meta$EQTL = 0
meta$EQTL[1:Ntop] = 1
meta.genes = meta$GENE[meta$EQTL == 1]

# Read in GWAS dataset
gwas = read.delim(paste0(annodir, '/GWAS/gwas_catalog_v1.0-downloaded_2015-11-30.tsv'), header = T, sep = '\t', stringsAsFactors = F)

# Get list of reported genes
reported.genes = unique(sort(unlist(strsplit(gwas$REPORTED.GENE.S., ' - |, '))))

# Read in OMIM dataset
omim.genes = read.table(paste0(annodir, '/OMIM/omim.genes.txt'), sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Read in OrphaNet dataset
orpha.genes.ensg = read.table(paste0(annodir, '/OrphaNet/orphanet.genes.txt'), sep = '\t', header = F, stringsAsFactors = F)[, 1]
orpha.genes = ensg.hgnc[orpha.genes.ensg, 1]
orpha.genes = orpha.genes[!is.na(orpha.genes)]

# Read in Cardio and Cancer genes 
cardio.gold.genes = read.table(paste0(annodir, '/GeneListsOther/cardio.genes.gold.standard.csv'), sep = ',', header = T, stringsAsFactors = F)[, 1]
cancer.gold.genes = read.table(paste0(annodir, '/GeneListsOther/cancer.genes.gold.standard.csv'), sep = ',', header = T, stringsAsFactors = F)[, 1]

# Read in ACMG data
acmg = read.table(paste0(annodir, '/ACMG/completed_acmg.csv'), sep = ',', header = T, stringsAsFactors = F)
acmg.genes = acmg$Gene

# Read in and define PTV-constrained genes (pLI > 0.9)
exac.threshold = 0.9
exac = read.table(paste0(exacdir, '/forweb_cleaned_exac_r03_march16_z_data_pLI.txt'), sep = '\t', stringsAsFactors = F, header = T)
exac.genes = exac$gene[exac$pLI > exac.threshold]

# Read in the ClinVar genes
clin = read.csv(paste0(annodir, '/GeneListsOther/gene_condition_source_id'), sep = '\t', stringsAsFactors = F)
clin.genes = clin$GeneSymbol

# Make list of gene lists to test
gene.lists = list(meta.genes,
	reported.genes, 
	omim.genes, 
	orpha.genes,
	clin.genes,  
	acmg.genes, 
	cardio.gold.genes, 
	cancer.gold.genes, 
	exac.genes)

# Calculate ORs for outlier and shared eQTL genes
medz.ors = gene.ors.calc(medz, 'Outlier', 'OUTLIER', gene.lists, gene.lists.names)
meta.ors = gene.ors.calc(meta, 'eGene', 'EQTL', gene.lists[-1], gene.lists.names[-1])

# Write out table of counts
write.table(medz.ors$table, paste0(dir, '/data/medz.disease.genes.table.txt'), sep = '\t', row.names = F, col.names = T, quote = F)
write.table(meta.ors$table, paste0(dir, '/data/metasoft.disease.genes.table.txt'), sep = '\t', row.names = F, col.names = T, quote = F)

# Combine the OR dataframes for easy plotting
full.ors = rbind(medz.ors$ors[-which(medz.ors$ors == 'eGene'), ], meta.ors$ors)
full.ors$TYPE = factor(full.ors$TYPE, levels = c('eGene', 'Outlier'))

# Plot the enrichments
# Remove GTEx eGenes
bad.lists = c(gene.lists.names[1], gene.lists.names[8])

medz.ors.plot = ggplot(data = medz.ors$ors[!(medz.ors$ors$SET %in% bad.lists[-length(bad.lists)]), ], aes(x = SET, y = OR, colour = SET)) + geom_abline(intercept = 1, slope = 0, colour = 'darkgrey') + 
	geom_pointrange(aes(x = SET, ymin = LOW, ymax = HIGH), size = 0.3) + theme_bw() + coord_flip() + ylab('Odds ratio') + xlab('') + 
	scale_colour_manual(values = colors) + guides(colour = F)
medz.ors.plot

alt.colors = c(colors[1], 'dodgerblue3')
names(alt.colors) = c('eGene', 'Outlier')

full.ors.plot.alt = ggplot(data = full.ors[!(full.ors$SET %in% bad.lists), ], aes(x = SET, y = OR, colour = TYPE)) +
    geom_abline(intercept = 1, slope = 0, colour = 'darkgrey') + 
    coord_flip() + geom_pointrange(aes(x = SET, ymin = LOW, ymax = HIGH), size = .3, position = position_dodge(width = 0.7)) +
    theme_bw() + ylab('Odds ratio') + xlab('') + 
    scale_colour_manual(values = alt.colors) +
    guides(colour = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position = c(.8, .9),
          legend.key = element_blank(),
          legend.text = element_text(size = 10),
          legend.background = element_blank())
full.ors.plot.alt

# Write out information for each disease list
write.gene.table('gwas.reported.outliers.txt', reported.genes)
write.gene.table('omim.outliers.txt', omim.genes)
write.gene.table('orpha.outliers.txt', orpha.genes)
write.gene.table('clinvar.outliers.txt', clin.genes)
write.gene.table('acmg.outliers.txt', acmg.genes)
write.gene.table('cardio.outliers.txt', cardio.gold.genes)
write.gene.table('cancer.outliers.txt', cancer.gold.genes)
write.gene.table('exac.lof.outliers.txt', exac.genes)

# Save workspace image
save.image(paste0(dir, '/data/figure4c.gene.list.enrichments.RData')




