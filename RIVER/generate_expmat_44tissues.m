dirNames.RAREVARDIR = '${RAREVARDIR}'

load(sprintf('%s/preprocessing/dataGTEx_052516.mat',dirNames.RAREVARDIR));
load(sprintf('%s/preprocessing/annoGTEx_052516.mat',dirNames.RAREVARDIR));
list_tissues = data.tisNames;

nGenes = length(anno.chr);
nInds = length(data.indNames.wgs);

target_data.exp = [];

matrix.gene_idx = repmat([1:nGenes]',1,nInds);
matrix.ind_idx = repmat([1:nInds],nGenes,1);
target_data.all = [reshape(matrix.gene_idx,nGenes*nInds,1) reshape(matrix.ind_idx,nGenes*nInds,1)];
for tis = 1:length(list_tissues)    
    % You might need to check if you can succefully upload expression data and convert missing values into NaN for matlab environments
    temp_exp = importdata(sprintf('%s/preprocessing/peer/%s.peer.ztrans.txt',dirNames.RAREVARDIR,list_tissues{tis}));
    if size(temp_exp.data,2) ~= size(temp_exp.textdata,2)-1,
        temp_exp.data = [temp_exp.data NaN(size(temp_exp.data,1),size(temp_exp.textdata,2)-size(temp_exp.data,2)-1)];
    end
    
    gene2ind = NaN(nGenes,nInds);    
    
    idx_wgs = NaN(2,length(data.indNames.wgs)); % new | old 
    for i = 1:length(data.indNames.wgs)
        idx_wgs(:,i) = [i; strmatch(data.indNames.wgs(i),temp_exp.textdata(1,2:end),'exact')];
    end
    
    idx_gene = NaN(length(anno.chr),2); % new | old
    for i = 1:length(anno.chr)
        temp_idx = strmatch(anno.gene_ids(i,2),temp_exp.textdata(2:end,1),'exact');
        if length(temp_idx) == 1,
            idx_gene(i,:) = [i temp_idx];
        end
    end
    idx_gene(sum(isnan(idx_gene),2) > 0,:) = [];
    
    gene2ind(idx_gene(:,1),idx_wgs(1,:)) = temp_exp.data(idx_gene(:,2),idx_wgs(2,:));
    save(sprintf('%s/RIVER/data/expression/%s.mat',dirNames.RAREVARDIR,list_tissues{tis}),'gene2ind');
    
    disp([' *** ' num2str(tis) 'th tissue data was generated and uploaded! *** ']);
    
    target_data.exp = [target_data.exp reshape(gene2ind,nGenes*nInds,1)]; 
end
value.median = nanmedian(target_data.exp,2);
value.nTissue = sum(~isnan(target_data.exp),2); clear gene2ind 
gene2ind.median = reshape(value.median,nGenes,nInds);
gene2ind.nTissue = reshape(value.nTissue,nGenes,nInds);
save(sprintf('%s/RIVER/data/expression/exp_median.mat',dirNames.RAREVARDIR),'gene2ind');
