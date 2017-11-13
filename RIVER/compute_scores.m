order = 1;
%coding = 1;
% =================================================================
% [step 2] Standardization % rescaling with imputed values (0)
% =================================================================

% Extract various features per gene per indiv
addpath(genpath('{RAREVARDIR}/RIVER/code'));
dirNames.RAREVARDIR = '${RAREVARDIR}';

load(sprintf('%s/reference/annoGTEx_052516.mat',dirNames.RAREVARDIR));
load(sprintf('%s/reference/dataGTEx_052516.mat',dirNames.RAREVARDIR));

process_nfeatures = 10;

% switch coding,
%     case 1,
%         option_rv = 'w_pc';
%     case 0,
%         option_rv = 'wo_pc';
% end

load(sprintf('%s/RIVER/data/score/indiv/GTEX-N7MS.1.score.gene.mat',dirNames.RAREVARDIR))
list_features = features.names;
% % list_features.cat_regions = zeros(length(features.names),1);
% % for i = 1:length(features.names)
% %     temp_names = strsplit(features.names{i},'_');
% %     switch temp_names{end},
% %         case 'tss3k',
% %             list_features.cat_regions(i) = 1;
% %         case 'tss5k',
% %             list_features.cat_regions(i) = 2;
% %         case 'tss10k',
% %             list_features.cat_regions(i) = 3;
% %         case 'tss20k',
% %             list_features.cat_regions(i) = 4;
% %         case 'tss50k',
% %             list_features.cat_regions(i) = 5;
% % %         case 'tss100k',
% % %             list_features.cat_regions(i) = 6;
% % %         case 'tss200k',
% % %             list_features.cat_regions(i) = 7;
% %     end
% % end
norm_info.mean = zeros(length(list_features),1);
norm_info.std = zeros(length(list_features),1);
save(sprintf('%s/RIVER/data/score/list_features_all.mat',dirNames.RAREVARDIR),'list_features','norm_info')

load(sprintf('%s/RIVER/data/score/list_features_all.mat',dirNames.RAREVARDIR))

nUniqFeatures = length(list_features); % tss10k only

label_region = {'10kb'};
% label_region = {'tss10k'};
count_score = 0;

for idx_score = 1+process_nfeatures*(order-1):process_nfeatures*(order)
    
    % tic
    count_score = count_score + 1;
    
    %% Raw features
    for r = 1:length(label_region)
        score(r).matrix = NaN(length(anno.chr),length(data.indNames.wgs));
    end
    n_gene = size(score(1).matrix,1);
    n_ind = size(score(1).matrix,2);
    
    for ii = 1:length(data.indNames.wgs)        
        load(sprintf('%s/%s.%s.mat',dirNames.score_src,data.indNames.wgs{ii},option_rv));        
        
        for r = 1:length(label_region)
            score(r).matrix(features.region.idx_genes,ii) = features.region.values(:,idx_score);
        end
        
        if rem(ii,10) == 0,
            disp([' *** ' num2str(ii) 'th indiv data is uploaded *** ']);
        end
    end
    
    for r = 1:length(label_region)
        gene2ind = score(r).matrix;        
        save(sprintf('%s/RIVER/data/score/feature/%s_raw.mat',dirNames.RIVERVARDIR,features.names{idx_score+nUniqFeatures*(r-1)}),'gene2ind');        
    end
    
    %% Scaled features with imputed values as 0
    % distance to TSS should be reversed before standardizing them
    
    temp_names = strsplit(features.names{idx_score},'_');
    for r = 1:length(label_region)
%         if strcmp(temp_names{1},'min') == 1, % min dist to TSS
%             switch label_region{r},
% %                 case 'tss3k',   ref_dist = 3e3+3e2;
% %                 case 'tss5k',   ref_dist = 5e3+5e2;
%                 case '10kb',  ref_dist = 1e4+1e3;
% %                 case 'tss20k',  ref_dist = 2e4+2e3;
% %                 case 'tss50k',  ref_dist = 5e4+5e3;
%                     %                 case 'tss100k', ref_dist = 1e5+1e4;
%                     %                 case 'tss200k', ref_dist = 2e5+2e4;
%             end
%             score(r).matrix = ref_dist - score(r).matrix;
%         end
        
        % standardization
        X = reshape(score(r).matrix,n_gene*n_ind,1);
        idx_nonan = find(isnan(X) == 0);
        x_sub = X(idx_nonan); x_sub = x_sub';
        [x_sub,m,st] = standardize(x_sub);
        X(idx_nonan) = x_sub';
        X = X - nanmin(X) + 5;
        
        norminfo = [m st];        
        
        score(r).matrix = reshape(X,n_gene,n_ind);
        score(r).matrix(isnan(score(r).matrix)) = 0;
        
        gene2ind = score(r).matrix;
        
        save(sprintf('%s/RIVER/data/score/feature/%s_scaled.mat',dirNames.RAREVARDIR,features.names{idx_score+nUniqFeatures*(r-1)}),'gene2ind','norminfo');        
    end
    disp([' === ' num2str(count_score) 'th feature: ' features.names{idx_score} ' === ']);
    %     toc
end
