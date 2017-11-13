%%%%% Assign directories into your own directories 
dirNames.RAREVARDIR = '${RAREVARDIR}';
%%%%%

% GTEx annotation file
fid = fopen(sprintf('%s/reference/gencode.v19.genes.v6p.patched_contigs.coding.lincRNA.gtf',dirNames.RAREVARDIR),'r');
data_anno = textscan(fid,'%d%d%d%s%s%s%s');
fclose(fid);
temp_anno.chr = double(data_anno{1});
temp_anno.pos = [double(data_anno{2}) double(data_anno{3})];
temp_anno.gene_ids = [data_anno{4} data_anno{5} data_anno{6} data_anno{7}];
% save annoGTEx_052516.mat anno

% GTEx data file
data.tisNames = importdata(sprintf('%s/reference/tissue_names.txt',dirNames.RAREVARDIR));
data.indNames.wgs = importdata(sprintf('%s/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt',dirNames.RAREVARDIR));
save(sprintf('%s/reference/dataGTEx.mat',dirNames.RAREVARDIR),'data')

% Expression Ensembl_ID file
ensembl_ids = importdata(sprintf('%s/reference/gene_ensembl_ids.txt',dirNames.RAREVARDIR));
target_ensembl_ids = intersect(ensembl_ids,temp_anno.gene_ids(:,2));

%fid = fopen('target_gene_ensembl_ids.txt','w');
%for i = 1:length(target_ensembl_ids)
%    fprintf(fid,'%s\n',target_ensembl_ids{i});
%end
%fclose(fid);

idx_target = zeros(length(target_ensembl_ids),1);
for i = 1:length(target_ensembl_ids)
    idx_target(i) = strmatch(target_ensembl_ids(i),temp_anno.gene_ids(:,2),'exact');
end
anno.chr = temp_anno.chr(idx_target);
anno.pos = temp_anno.pos(idx_target,:);
anno.gene_ids = temp_anno.gene_ids(idx_target,:);
save(sprintf('%s/reference/annoGTEx.mat',dirNames.RAREVARDIR),'anno')                                                                                                                                                  






