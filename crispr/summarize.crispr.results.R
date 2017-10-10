#!/usr/bin/R

rm(list = ls())

dir = Sys.getenv('RAREVARDIR')

# Load packages
require(ggplot2)
require(cowplot)
require(ggbeeswarm)
require(reshape2)
require(plyr)
require(dplyr)
require(broom)
require(gridExtra)
require(Rsamtools)
require(stringr)
require(foreach)
require(doMC)

registerDoMC(cores = 10)

#----------- FUNCTIONS

sum.agg <- function(x){
	out = sum(x, na.rm = T)
	return(out)
}

#----------- MAIN

# Global constants
types = c('Outlier', 'Control')

dna.type.colors = c('darkgrey', 'dodgerblue4')
dna.type = c('gDNA', 'cDNA')
names(dna.type.colors) = dna.type

fsize = 11
gtex.theme = theme_bw() + theme(axis.text = element_text(size = fsize),
	axis.title = element_text(size = fsize), 
	legend.text = element_text(size = fsize), 
	legend.title = element_text(size = fsize), 
	strip.text = element_text(size = fsize + 1), 
	strip.background = element_blank(), 
	plot.title = element_text(hjust = .5, size = fsize + 1), 
	legend.background = element_blank(),
	legend.key = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank())

count.scale = 	scale_size_discrete(labels = c(expression(10^2), expression(10^3), expression(10^4), 
	expression(10^5), expression(10^6)), name = 'Read count')

# Define the experimental design matrix
vars = c('11:64956213:C:G', '7:66459273:T:A', '12:4766944:C:T', '7:102944937:G:A', '19:13885293:T:A', '3:120321157:A:G', 
	'11:64972286:G:C', '7:66459256:T:C', '12:4766925:G:T', '7:102948074:A:G', '19:13885309:C:T')
gdna.out.pos = c(183, 209, 168, 212, 172, 161)
cdna.out.pos = c(158, 165, 149, 198, 37, 219)
gdna.cont.pos = c(164, 192, 149, 179, 167)
cdna.cont.pos = c(295, 182, 130, 230, 53)
exp.design = data.frame(Sample = rep(1:66, 2), 
	Rep = rep(c(rep(rep(1:3, each = 6), 2), rep(1:3, 10)), 2),  
	Type = factor(rep(c(rep(types[1], 36), rep(types[2], 30)), 2), levels = types), 
	DNA = factor(rep(c(rep(dna.type, each = 18), rep(dna.type, each = 15)), 2), levels = dna.type), 
	RefAlt = factor(rep(c('Ref', 'Alt'), each = 66), levels = c('Ref', 'Alt')), 
	Pos = rep(c(rep(gdna.out.pos, 3),
		rep(cdna.out.pos, 3), 
		rep(gdna.cont.pos, each = 3), 
		rep(cdna.cont.pos, each = 3)), 2), 
	VarId = factor(rep(c(rep(vars[1:6], 6),
		rep(rep(vars[7:11], each = 3), 2)), 2), levels = vars))
exp.design = exp.design %>% mutate(ID = paste(DNA, VarId, RefAlt, Pos, sep = '_'))

write.table(exp.design, paste0(dir, '/data/CRISPR/crispr.design.matrix.txt'),
            row.names = F, col.names = T, quote = F, sep = "\t")

# Make pileups for each sequencing run 
# Filter for loci of interest
bam.files = list.files(paste0(dir, '/data/CRISPR/bams', pattern = '_n0.bam', full.names = T)
pileups = foreach(i=1:length(bam.files), .combine = rbind) %dopar% {
	print(i)
	sample.id = as.numeric(gsub('_n0.bam', '', gsub('.*/SBM_Samp', '', bam.files[i])))
	locus.id = exp.design %>% filter(Sample == sample.id)
	sample.pileup = pileup(BamFile(file = bam.files[i], index = gsub('_n0.bam', '_n0.bai', bam.files[i])),
		pileupParam = PileupParam(max_depth = 1e9, min_mapq = 20, min_base_quality = 0, min_nucleotide_depth = 0, distinguish_strands = F))
	sample.pileup = sample.pileup %>% mutate(Sample = sample.id, Rep = locus.id$Rep[1], Type = locus.id$Type[1], VarId = locus.id$VarId[1], 
		ID = paste(seqnames, pos, sep = '_'), RefAlt = gsub('.*_', '', seqnames), DNA = gsub('_.*', '', seqnames), 
		Allele = str_split(VarId, ':')[[1]][ifelse(RefAlt == 'Ref', 3, 4)]) %>%
		filter(ID %in% locus.id$ID) %>%
		select(Sample, Rep, Type, DNA, RefAlt, Allele, pos, VarId, nucleotide, count)
	names(sample.pileup)[c(7, 9)] = c('Pos', 'Allele.Pileup')
	sample.pileup
}
pileups = pileups[order(pileups$Sample), ]
pileups = merge(exp.design, pileups, all.x = T)
refAlt = pileups %>% dcast(., Sample + Rep + Type + DNA + VarId ~ RefAlt, fill = 0, fun.aggregate = sum.agg, value.var = 'count') %>%
	mutate(AltFreq = Alt / (Ref + Alt), FC = Alt / Ref, Count = Ref + Alt)

# Read in Picard summary stats on total reads, aligned reads, adapter contamination, etc.
# Merge with the allele count data
picard.files = list.files(paste0(dir, '/data/CRISPR/bams'), pattern = 'summ.stats', full.names = T)
picard = foreach(i=1:length(picard.files), .combine = rbind) %dopar%{
	print(i)
	sample.id = as.numeric(gsub('_n0.*', '', gsub('.*/SBM_Samp', '', picard.files[i])))
	sample.info = exp.design %>% filter(Sample == sample.id) %>% filter(RefAlt == 'Ref')
	stats = read.table(picard.files[i], sep = '\t', header = T, stringsAsFactors = F) %>%
		select(TOTAL_READS, PF_READS, PCT_PF_READS, PF_READS_ALIGNED, PCT_PF_READS_ALIGNED, PF_HQ_ALIGNED_READS, MEAN_READ_LENGTH) %>%
		mutate(Sample = sample.info$Sample[1], Rep = sample.info$Rep[1], Type = sample.info$Type[1], DNA = sample.info$DNA[1], VarId = sample.info$VarId[1])
	stats
}
refAlt = merge(refAlt, picard, all = T)
refAlt = refAlt[order(refAlt$Sample), ]

# Write out Alt/Ref data for sharing
write.table(refAlt, paste0(dir, '/data/CRISPR/crispr.coding.control.sites.ref.alt.gdna.cdna.txt'), sep = '\t', col.names = T, row.names = F, quote = F)

# Remove the first and last outlier loci and the first control locus because they have low read depth
refAlt = refAlt %>% filter(!(VarId %in% c('3:120321157:A:G', '11:64956213:C:G', '11:64972286:G:C')))

# Make a factor corresponding to the read count to use for the point size
refAlt = refAlt %>% mutate(Size = cut(Count, breaks = 10^c(2:7), labels = 2:6))

# Make a plot showing the experimental design
exp.design = exp.design %>% filter(RefAlt == 'Ref' & VarId != '3:120321157:A:G') %>%
	select(Rep, Type, DNA, VarId) %>% mutate(Rep = ifelse(DNA == 'gDNA', Rep, Rep + 3))

ed.outlier = ggplot(data = filter(exp.design, Type == 'Outlier'), aes(x = Rep, y = VarId)) +
	gtex.theme + geom_tile(aes(fill = DNA), colour = 'white', size = 2) +
	scale_fill_manual(values = dna.type.colors) + guides(fill = F) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), panel.border = element_blank(), plot.title = element_text(hjust = .5), 
		axis.ticks = element_blank()) +
	xlab('Sample') + ylab('') + ggtitle('Outlier') + scale_x_continuous(breaks = 1:6)

ed.control = ggplot(data = filter(exp.design, Type == 'Control'), aes(x = Rep, y = VarId)) +
	gtex.theme + geom_tile(aes(fill = DNA), colour = 'white', size = 2) +
	scale_fill_manual(values = dna.type.colors) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), panel.border = element_blank(), plot.title = element_text(hjust = .5), 
		axis.ticks = element_blank()) +
	xlab('Sample') + ylab('') + ggtitle('Control') + scale_x_continuous(breaks = 1:6)

ed.plot = ggdraw() + 
  draw_plot(ed.outlier, 0, 0, .45, 1) +
  draw_plot(ed.control, .45, 0, .55, 1) 
ed.plot

# Test the difference in gDNA and cDNA alt allele proportions
# Use a simple t-test
t.tests = refAlt %>% group_by(VarId, Type) %>%
	do(tidy(t.test(AltFreq ~ DNA, data = .)))
t.tests$FDR = p.adjust(t.tests$p.value, method = 'bonferroni')
t.tests$stars = cut(t.tests$FDR, breaks = c(1, 0.05, 1e-2, 1e-3, 1e-4, 1e-5), labels = rev(c('', '.', '*', '**', '***')))
t.tests$Y = unlist(lapply(t.tests$VarId, function(x){max(refAlt$AltFreq[refAlt$VarId == x]) + .05}))

# Make plot of alt allele frequency for gDNA and cDNA for each locus
locus.ids = data.frame(VarId = vars[c(2:5, 8:11)], LocusId = rep(1:4, 2))
refAlt = merge(refAlt, locus.ids, all.x = T)
t.tests = merge(t.tests, locus.ids, all.x = T)

alt.freq.plot = ggplot(data = refAlt, aes(x = LocusId, y = AltFreq)) +
	theme_bw() + gtex.theme + geom_quasirandom(aes(colour = DNA, size = Size), width = .4) +
	count.scale +
	xlab('Locus ID') + ylab('Alternate allele proportion') +
	scale_colour_manual(values = dna.type.colors, name = '(g/c)DNA') +
	geom_text(data = t.tests, aes(x = LocusId, y = Y, label = stars), size = 8) +
	guides(colour = guide_legend(override.aes = list(size = 5.5), order = 1), size = guide_legend(order = 2)) +
	ylim(c(0, .45)) + facet_grid(. ~ Type)
alt.freq.plot

#---- Make plots of total read count, aligned reads, etc. 
# Total read count
trc.plot = ggplot(data = refAlt, aes(x = LocusId, y = TOTAL_READS)) +
	theme_bw() + gtex.theme + geom_quasirandom(aes(colour = DNA), alpha = .5, width = .4, size = 5.5) +
	scale_y_log10(limits = c(1e6, 1e8)) +
	xlab('Locus ID') + ylab('Total read counts') +
	scale_colour_manual(values = dna.type.colors, name = '(g/c)DNA') +
	facet_grid(. ~ Type)
trc.plot

# Calculate median read depth per sample
median(refAlt$TOTAL_READS)
# 4710266

# Percent of total reads that are high quality aligned reads
hqa.plot = ggplot(data = refAlt, aes(x = LocusId, y = PF_HQ_ALIGNED_READS / TOTAL_READS)) +
	theme_bw() + gtex.theme + geom_quasirandom(aes(colour = DNA), alpha = .5, width = .4, size = 5.5) +
	xlab('Locus ID') + ylab('Percent high quality aligned reads (MAPQ > 20)') +
	scale_colour_manual(values = dna.type.colors, name = '(g/c)DNA') +
	ylim(c(0, .3)) +
	facet_grid(. ~ Type)
hqa.plot

# Calculate median high quality aligned reads per sample
median(refAlt$PF_HQ_ALIGNED_READS / refAlt$TOTAL_READS)
# 0.08557442

# Calculate it for the cDNA for 7:66459256:T:C
refAlt %>% filter(DNA == 'cDNA' & VarId == '7:66459256:T:C') %>% summarize(median(PF_HQ_ALIGNED_READS / TOTAL_READS))
#  median(PF_HQ_ALIGNED_READS/TOTAL_READS)
#1                            0.0001210156

# number of high quality aligned reads
nhq.plot = ggplot(data = refAlt, aes(x = LocusId, y = PF_HQ_ALIGNED_READS)) +
	theme_bw() + gtex.theme + geom_quasirandom(aes(colour = DNA), alpha = .5, width = .4, size = 5.5) +
	scale_y_log10(limits = c(1, 1e7)) +
	xlab('Locus ID') + ylab('Total high quality aligned reads (MAPQ > 20)') +
	scale_colour_manual(values = dna.type.colors, name = '(g/c)DNA') +
	facet_grid(. ~ Type)
nhq.plot

# Percent of high quality reads that align to the locus of interest
phq.plot = ggplot(data = refAlt, aes(x = LocusId, y = Count / PF_HQ_ALIGNED_READS)) +
	theme_bw() + gtex.theme + geom_quasirandom(aes(colour = DNA), alpha = .5, width = .4, size = 5.5) +
	xlab('Locus ID') + ylab('Percent high quality aligned reads\nmapping to the locus of interest') +
	scale_colour_manual(values = dna.type.colors, name = '(g/c)DNA') +
	guides(colour = F) + ylim(c(.7, 1)) +
	facet_grid(. ~ Type)
phq.plot

# Save workspace
save.image(paste0(dir, '/data/summarize.crispr.results.RData'))

