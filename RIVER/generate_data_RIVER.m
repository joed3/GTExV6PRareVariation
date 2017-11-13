mask_missing = 0;
per_missing = 0.5;
median_thrd.gene = 2;
median_thrd.outlier = 1.5;
label_median_thrd.gene = '2';
label_median_thrd.outlier = '1p5';

dirNames.RAREVARDIR = '${RAREVARDIR}'

% annotation files
load(sprintf('%s/reference/dataGTEx_052516.mat',dirNames.RAREVARDIR)); data_info.nIndivs = length(data.indNames.wgs);
load(sprintf('%s/reference/annoGTEx_052516.mat',dirNames.RAREVARDIR)); data_info.nGenes = length(anno.chr);

num_indivs = length(data.indNames.wgs);
num_gene = length(anno.chr);

list_tissues = data.tisNames;
ntissues = length(list_tissues);

load(sprintf('%s/RIVER/data/score/list_features_all.mat',dirNames.RAREVARDIR));
score_list = list_features;%(find(list_features.cat_regions == target_region));

for coding = 1%[0 1]; % 0: rare variants only in non-coding regions
    switch coding,
        case 0, option_rv = 'wo_pc';
        case 1, option_rv = 'w_pc';
    end
    
    for mask_missing = 0           
        disp([' *** [missing_mask,coding] = [' num2str(mask_missing) ',' num2str(coding) '] *** ']);        
        
        % find maximum # of indivs having rare variants (10K from TSS)
        idx_tss10k = [1:length(list_features)];%find(list_features.cat_regions == target_region);
        idx_target_feature = 0;
        idx_examples_rare = 0;
        target_nonan = 0;
        for nScore = 1:length(idx_tss10k)
            load(sprintf('%s/RIVER/data/score/feature/%s_scaled.mat',dirNames.RAREVARDIR,list_features{idx_tss10k(nScore)}));
            score = gene2ind';
            temp_bg = 0; % background value (always 0)
            
            score(score == temp_bg) = NaN;
            n_nonan = sum(sum(~isnan(score)));
            if n_nonan > target_nonan;
                target_nonan = n_nonan;
                % index of g-feature having highest # of indivs with rare variants
                idx_target_feature = idx_tss10k(nScore);
            end
        end
        
        %% genomic mask based on a genomic feature having max. # of rare variants
        % values: feature values with a rare variant
        % NaN   : either missing feature or no rare variant
        load(sprintf('%s/RIVER/data/score/feature/%s_scaled.mat',dirNames.RAREVARDIR,list_features{idx_target_feature}));
        mask_data.score = gene2ind';
        temp_bg = 0; % background value (always 0)
        mask_data.score(mask_data.score == temp_bg) = NaN;
        disp([' ++++++ Overall genomic mask was generated ... ++++++ ']);
        
        %% expression masks
        mask_data.exp = [];
        for tis = 1:length(list_tissues)
            load(sprintf('%s/RIVER/data/expression/%s.mat',dirNames.RAREVARDIR,list_tissues{tis}));
            %     load(sprintf('%s/%s.wgs.bin.mat',dirNames.e,CorrMatrix.names{tis}));
            ind2gene = gene2ind';
            mask_data.exp = [mask_data.exp reshape(ind2gene,data_info.nGenes*data_info.nIndivs,1)];
        end
        % % [outliers] # of outliers across corresponding tissues per instance
        % % values: # of outliers across tissues
        % mask_data.outliers = nansum(mask_data.exp,2);
        % mask_data.outliers = reshape(mask_data.outliers,data_info.nIndivs,data_info.nGenes);
        
        n_tissue_thrd = 5;
        load(sprintf('%s/RIVER/data/expression/exp_median.mat',dirNames.RAREVARDIR));
        median_data.exp = gene2ind; clear gene2ind
        
        ind2gene.median = median_data.exp.median';
        ind2gene.nTissue = median_data.exp.nTissue'; clear median_data;
        ind2gene.median(ind2gene.nTissue < n_tissue_thrd) = NaN;    % filter out any examples having expressions < n_tissue_thrd
        
        % processing expression data
        temp_exp.org = abs(ind2gene.median);
        temp_exp.bin_sin = NaN(size(temp_exp.org));
        temp_exp.bin_mul = NaN(size(temp_exp.org));
        
        count_outlier = 0;
        n_outliers.median = NaN(1,size(temp_exp.org,2));
        for g = 1:size(temp_exp.org,2)
            %     temp_vec = ind2gene.median(:,g);
            idx_nonan = find(isnan(temp_exp.org(:,g)) == 0);
            temp_vec = temp_exp.org(idx_nonan,g);
            temp_vec1 = tiedrank(temp_vec);
            
            temp_vec_sin = zeros(length(idx_nonan),1);
            idx_max_sin = find((temp_vec >= median_thrd.gene) & (temp_vec1 == max(temp_vec1)));
            if length(idx_max_sin) > 0,
                count_outlier = count_outlier + 1;
                temp_vec_sin(idx_max_sin) = 1;
            end
            temp_exp.bin_sin(idx_nonan,g) = temp_vec_sin;
            
            temp_vec_mul = zeros(length(idx_nonan),1);
            idx_max_mul = find(temp_vec >= median_thrd.outlier);
            
            % outlier gene: |median(z)| >= 2, outlier calling: multiple
            if length(idx_max_sin) > 0,
                n_outliers.median(g) = length(idx_max_mul);
                temp_vec_mul(idx_max_mul) = 1;
            end
            temp_exp.bin_mul(idx_nonan,g) = temp_vec_mul;
            %     temp_vec = max(temp_vec) - temp_vec + 1;
            %     ind2gene.median(idx_nonan,g) = temp_vec2;
        end
%         mask_data.outliers = temp_exp.bin_sin;
        mask_data.outliers = temp_exp.bin_mul; % simple threshold based outlier calling (|median(z)| > K)
        
        disp([' ... median expression data was processed']);
        disp(['     # outliers = ' num2str(count_outlier)]);
        % [bin_s, bin_m, cont] = Get_values_stouffer(ntissues,list_tissues,data_info,dirNames);
        % mask_data.outliers = bin_s;
        
        % find target genes having at least one outlier indiv. in any of tissues due to rare variants
        idx_noout = find(nansum(mask_data.outliers,1) == 0);
        mask_data.outliers(:,idx_noout) = NaN;
        
        class_labels.ind2gene.sin = temp_exp.bin_sin;
        class_labels.ind2gene.mul = temp_exp.bin_mul;
        % class_labels.ind2gene.soft = cont;  clear bin_s bin_m cont
        
        % [missing] # of missing values across corresponding tissues per instance
        % values: # of missing values across tissues
        % NaN: percent missing values > 50 %
        mask_data.missing = sum(isnan(mask_data.exp),2);
        mask_data.missing = reshape(mask_data.missing,data_info.nIndivs,data_info.nGenes);
        mask_data.missing(mask_data.missing >= ntissues*per_missing) = NaN; % percent missingness < 50 %
        
        %% combined mask
        switch mask_missing,
            case 1, % filtering instances with % missingness
                mask_data.combined = mask_data.score.*mask_data.outliers.*mask_data.missing;
            case 0, % filtering instances without % missingness
                mask_data.combined = mask_data.score.*mask_data.outliers;
        end
        
        % list of indices for target ind x gene instances
        % mask_target = ~isnan(mask_data.combined).*~isnan(mask_data.sel_genes); clear mask_data
        mask_target = ~isnan(mask_data.combined); %clear mask_data
        [ind, gene] = find(mask_target == 1);
        % combined.idx_ind_gene = [ind gene]; clear ind gene
        combined.names.ind_gene = [data.indNames.wgs(ind) anno.gene_ids(gene,2)]; clear ind gene        
        
        N = size(combined.names.ind_gene,1); % # of indiv-gene instances
        P = length(score_list);             % # of genomic features
        T = ntissues;                       % size(combined.values.e,2);            % # of tissues
        num_nodes = T + 1;                  % T tissue nodes + 1 class node
        
        %% generate genomic feature data
        combined.names.g = score_list;
        combined.values.g = zeros(N,P);
        for nScore = 1:length(score_list)
	    load(sprintf('%s/RIVER/data/score/feature/%s_scaled.mat',dirNames.RAREVARDIR,score_list{nScore}));
            score = gene2ind';
            combined.values.g(:,nScore) = score(mask_target == 1);
        end
        
        temp_values = [];
        for tis = 1:length(list_tissues)
            load(sprintf('%s/RIVER/data/expression/%s.mat',dirNames.RAREVARDIR,list_tissues{tis}));
            ind2gene = gene2ind';
            temp_values = [temp_values ind2gene(mask_target == 1)];
        end
        
        %% generate expression data
        temp_bnet.names.e = list_tissues;
        temp_bnet.values.e_cont = temp_values;
%         temp_bnet.values.e_disc = temp_values;
%         temp_bnet.values.e_disc(abs(temp_bnet.values.e_disc) < 2) = 0;
%         temp_bnet.values.e_disc(abs(temp_bnet.values.e_disc) >= 2) = 1; % outlier
        
%         combined.values.e_disc = temp_bnet.values.e_disc;
        combined.values.e_cont = temp_bnet.values.e_cont;
        combined.names.e = temp_bnet.names.e;
        combined.values.e_median_single = class_labels.ind2gene.sin(mask_target == 1);      % binary & single outlier
        combined.values.e_median_multiple = class_labels.ind2gene.mul(mask_target == 1);      % binary & multiple outlier
        clear class_labels
        
        disp([' # samples = ' num2str(size(combined.values.g,1)) ', # outliers = ' num2str(sum(combined.values.e_median_multiple))]);
        
        switch mask_missing,
            case 1, % with % missingness filter
                % genomic features 
                fid = fopen(sprintf('%s/RIVER/data/score/genomic_features_%dperMissing.txt',dirNames.RAREVARDIR,per_missing*100),'w');
                fprintf(fid,'indiv\tgene\t');
                for c = 1:length(combined.names.g)
                    fprintf(fid,'%s\t',combined.names.g{c});
                end
                fprintf(fid,'\n');
                for r = 1:size(combined.names.ind_gene,1)
                    fprintf(fid,'%s\t%s\t',combined.names.ind_gene{r,1},combined.names.ind_gene{r,2});
                    for c = 1:size(combined.values.g,2)
                        fprintf(fid,'%f\t',combined.values.g(r,c));
                    end
                    fprintf(fid,'\n');
                end
                fclose(fid);
                
                % zscores 
                fid = fopen(sprintf('%s/RIVER/data/score/zscores_%dperMissing.txt',dirNames.RAREVARDIR,per_missing*100),'w');
                fprintf(fid,'indiv\tgene\t');
                for c = 1:length(combined.names.e)
                    fprintf(fid,'%s\t',combined.names.e{c});
                end
                fprintf(fid,'\n');
                for r = 1:size(combined.names.ind_gene,1)
                    fprintf(fid,'%s\t%s\t',combined.names.ind_gene{r,1},combined.names.ind_gene{r,2});
                    for c = 1:size(combined.values.e_cont,2)
                        fprintf(fid,'%f\t',combined.values.e_cont(r,c));
                    end
                    fprintf(fid,'\n');
                end
                fclose(fid);
            case 0, % without % missingness filter
                % genomic features 
                fid = fopen(sprintf('%s/RIVER/data/score/genomic_features.txt',dirNames.RAREVARDIR),'w');
                fprintf(fid,'indiv\tgene\t');
                for c = 1:length(combined.names.g)
                    fprintf(fid,'%s\t',combined.names.g{c});
                end
                fprintf(fid,'\n');
                for r = 1:size(combined.names.ind_gene,1)
                    fprintf(fid,'%s\t%s\t',combined.names.ind_gene{r,1},combined.names.ind_gene{r,2});
                    for c = 1:size(combined.values.g,2)
                        fprintf(fid,'%f\t',combined.values.g(r,c));
                    end
                    fprintf(fid,'\n');
                end
                fclose(fid);
                
                % zscores
		fid = fopen(sprintf('%s/RIVER/data/score/zscores.txt',dirNames.RAREVARDIR),'w');
                fprintf(fid,'indiv\tgene\t');
                for c = 1:length(combined.names.e)
                    fprintf(fid,'%s\t',combined.names.e{c});
                end
                fprintf(fid,'\n');
                for r = 1:size(combined.names.ind_gene,1)
                    fprintf(fid,'%s\t%s\t',combined.names.ind_gene{r,1},combined.names.ind_gene{r,2});
                    for c = 1:size(combined.values.e_cont,2)
                        fprintf(fid,'%f\t',combined.values.e_cont(r,c));
                    end
                    fprintf(fid,'\n');
                end
                fclose(fid);
        end
    end
end






