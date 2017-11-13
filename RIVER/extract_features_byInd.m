order = 1;
%coding = 1;

% coding = 1; % with protein coding, coding = 0; % without protein coding
% =================================================================
% [step 1] Extract various features per gene and indiv
% 1. remove any rare variants within exon regions.
% 2. extract all the features
% =================================================================

%%%%% Assign directories into your own directories                                                                                                                                                          
dirNames.RAREVARDIR = '${RAREVARDIR}';
%%%%%

load(sprintf('%s/reference/annoGTEx_052516.mat',dirNames.RAREVARDIR));
load(sprintf('%s/reference/dataGTEx_052516.mat',dirNames.RAREVARDIR));
list_indivs = data.indNames.wgs;               % list of indivs (116, compact ids)

%% You need to use numeric chromosome numbers here by removing "chr" in this file 
system('sed "s/chr//" ${RAREVARDIR}/reference/gencode.v19.annotation_coding.lincRNA_padded.bed > gencode.v19.annotation_coding.lincRNA_padded.chr_num.bed');
list_pcoding = importdata(sprintf('%s/reference/gencode.v19.annotation_coding.lincRNA_padded.bed',dirNames.RAREVARDIR));

process_nind = 1;

%switch coding,
%    case 1,
%        option_rv = 'w_pc';
%    case 0,
%        option_rv = 'wo_pc';
%end

% variant annotation (12 types)
% mapping.var_anno = {'StopChange'; 'SpliceVariant'; 'mRNAChange'; 'Nonsynonymous'; ...
%     'Synonymous'; 'TranscriptModifier'; 'UTR'; 'Intron'; 'Downstream'; ...
%     'Upstream'; 'Regulatory'; 'Intergenic'};
% mapping.var_anno = flipud(mapping.var_anno);

mapping.var_anno = {'stop_gained','stop_lost','splice_acceptor_variant','splice_donor_variant',...
    'transcript_ablation','frameshift_variant','start_lost','transcript_amplification',...
    'inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant',...
    'splice_region_variant','incomplete_terminal_codon_variant','stop_retained_variant',...
    'synonymous_variant','coding_sequence_variant','mature_miRNA_variant','non_coding_transcript_exon_variant',...
    'NMD_transcript_variant','non_coding_transcript_variant','5_prime_UTR_variant',...
    '3_prime_UTR_variant','intron_variant','downstream_gene_variant','upstream_gene_variant',...
    'TFBS_ablation','TFBS_amplification','TF_binding_site_variant','regulatory_region_ablation',...
    'regulatory_region_amplification','regulatory_region_variant','feature_elongation',...
    'feature_truncation','intergenic_variant','undefined'};
mapping.var_features = {'StopGained','StopLost','SpliceAcceptorVariant','SpliceDonorVariant',...
    'TranscriptAblation','FrameshiftVariant','StartLost','TranscriptAmplification',...
    'InframeInsertion','InframeDeletion','MissenseVariant','ProteinAlteringVariant',...
    'SpliceRegionVariant','IncompleteTerminalCodonVariant','StopRetainedVariant',...
    'SynonymousVariant','CodingSequenceVariant','MatureMiRNAVariant','NonCodingTranscriptExonVariant',...
    'NMDTranscriptVariant','NonCodingTranscriptVariant','5PrimeUTRVariant',...
    '3PrimeUTRVariant','IntronVariant','DownstreamGeneVariant','UpstreamGeneVariant',...
    'TFBSAblation','TFBSAmplification','TFBindingSiteVariant','RegulatoryRegionAblation',...
    'RegulatoryRegionAmplification','RegulatoryRegionVariant','FeatureElongation',...
    'FeatureTruncation','IntergenicVariant','Undefined'};

% segway (12 types)
mapping.segway = {'CTCF'; 'Dead'; 'Enh'; 'FAIRE'; 'GeneEnd'; 'GeneMid'; ...
    'GeneStr'; 'H3K9me1'; 'Low'; 'Reprd'; 'TF'; 'TSS'};

% chromHMM (10 types)
mapping.chromHMM = {'ActProm'; 'WeakProm'; 'PoisedProm'; 'StrongEnh'; 'WeakEnh'; ...
    'Insulator'; 'Txn'; 'WeakTxn'; 'Reprd'; 'Low'};

% For either anno or segway or chromHMM, there are two types of features, nvar and
% ratio (site)

% score_labels.operations = {'max','min','nvar','ratio'};
score_labels.operations = {'max','min','nvar'};
% score_labels.regions = {'tss3k','tss5k','tss10k','tss20k','tss50k','tss100k','tss200k'};
score_labels.regions = {'10kb'};

for ii = 1+process_nind*(order-1):process_nind*order%length(list_indivs)
    id_ind = list_indivs{ii};
    
    features.names = [];            % list of feature names (column)
    
    tic
    count_features = 0;
    for r = 1:length(score_labels.regions) %
        temp_data = importdata(sprintf('%s/RIVER/data/score/indiv/%s.%d.score.nuc.txt',dirNames.RAREVARDIR,id_ind,ii));
        %         temp_data.data = [cellfun(@str2double,temp_data.textdata(2:end,2:4)) NaN(size(temp_data.data,1),1) temp_data.data];
        
        temp_names = temp_data.textdata(1,3:end); % except for Ensembl ID and anno
        
        idx_info.chr = strmatch('Chrom',temp_names,'exact');
        idx_info.pos = strmatch('Pos',temp_names,'exact');
        idx_info.nvar = strmatch('nvar',temp_names,'exact');
        %         idx_info.var_anno = strmatch('anno',temp_names,'exact');
        
        target_names = temp_names(idx_info.nvar+1:end);
        
        % imputing missing values with 0 from encode
        count_imp_features = 0;
        for i = 1:length(target_names)
            idx_enc = strmatch('Enc',target_names(i));
            idx_TFBS = strmatch('TFBS',target_names(i));
            if length(idx_enc) + length(idx_TFBS) > 0,
                idx_target_feature = strmatch(target_names(i),temp_data.textdata(1,:),'exact')-2;
                temp_data.data(isnan(temp_data.data(:,idx_target_feature)),idx_target_feature) = 0;
                count_imp_features = count_imp_features + 1;
            end
        end
        
        features.region(r).values = [];             % list of features values with rare variants
        features.region(r).idx_genes = [];          % list of indexes for genes having rare variants based on annotation matrix
        features.region(r).pos_all = [];            % list of rare variant positions
        features.region(r).pos_oper = [];           % selected positions in each feature based on corresponding operations
        features.region(r).pos_fnames = [];
        features.region(r).so = [];                 % sequence ontology
        
        count_gene = 0;
        count_coding = 0;
        count_noncoding = 0;
        factor = 0;
        for g = 1:length(anno.chr) %%
            temp_idx.gene = strmatch(anno.gene_ids(g,2),temp_data.textdata(2:end,1),'exact');
            %             temp_idx.gene = find(temp_data.data(:,idx_info.idx_all) == anno.idx_all(g));
            if length(temp_idx.gene) > 0, % Is this gene cosidered?
                count_gene = count_gene + 1;
                
                features.region(r).idx_genes = [features.region(r).idx_genes; g];   % index from annoGTEx
                features.region(r).pos_all = [features.region(r).pos_all; ...
                    {temp_data.data(temp_idx.gene,[idx_info.chr idx_info.pos])}];   % chr, position
                features.region(r).so = [features.region(r).so; ...
                    {temp_data.textdata(temp_idx.gene+1,2)}];     % so
                
                target_data.features = temp_data.data(temp_idx.gene,idx_info.nvar+1:end);    % gene-specific feature matrix
                target_data.nvar = temp_data.data(temp_idx.gene,idx_info.nvar);
                target_data.pos = temp_data.data(temp_idx.gene,idx_info.pos);
                target_data.anno = temp_data.textdata(temp_idx.gene+1,2);
                
                temp_features = [];
                temp_pos = [];
                
                % annotation
                var_mat = zeros(length(temp_idx.gene),length(mapping.var_anno));
                for var1 = 1:length(target_data.anno)
                    temp_var_anno = unique(regexp(target_data.anno{var1},'&','split'));
                    for var2 = 1:length(temp_var_anno)
                        var_mat(var1,strmatch(temp_var_anno{var2},mapping.var_anno,'exact')) = ...
                            var_mat(var1,strmatch(temp_var_anno{var2},mapping.var_anno,'exact')) + target_data.nvar(var1);
                    end
                end
                var_mat = sum(var_mat,1);
                
                for so = 1:length(mapping.var_anno)
                    if count_gene == 1,
%                         count_features = count_features + 2;
                        count_features = count_features + 1;
                        features.names = [features.names; {[score_labels.operations{3} '_so:' ...
                            mapping.var_features{so} '_' score_labels.regions{r}]}];
%                             {[score_labels.operations{4} '_so:' ...
%                             mapping.var_features{so} '_' score_labels.regions{r}]}];
                    end
                    
%                     temp_features = [temp_features var_mat(so) var_mat(so)/sum(var_mat)];
                    temp_features = [temp_features var_mat(so)];
                end
                
                for s = 1:length(target_names) %%%
                    switch target_names{s},
                        case 'Segway',
                            ln = sum(~isnan(target_data.features(:,s)));
                            for seg = 1:size(mapping.segway,1) %%%%
                                if count_gene == 1,
%                                     count_features = count_features + 2;
                                    count_features = count_features + 1;
                                    features.names = [features.names; {[score_labels.operations{3} '_seg:' ...
                                        mapping.segway{seg} '_' score_labels.regions{r}]}];
%                                         {[score_labels.operations{4} '_seg:' ...
%                                         mapping.segway{seg} '_' score_labels.regions{r}]}];
                                end
                                
                                % # of rare variants for corresponding chromatin states
                                temp_idx.state = find(target_data.features(:,s) == seg);
                                if ln == 0, % No chromatin states
%                                     temp_features = [temp_features NaN NaN];
                                    temp_features = [temp_features NaN];
                                elseif (ln > 0) && (length(temp_idx.state) > 0),
                                    temp_features = [temp_features sum(target_data.nvar(temp_idx.state))]; ...
%                                         length(temp_idx.state)/ln];
                                elseif (ln > 0) && (length(temp_idx.state) == 0),
%                                     temp_features = [temp_features 0 0];
                                    temp_features = [temp_features 0];
                                end
                            end %%%%
                        case 'chromHMM',
                            ln = sum(~isnan(target_data.features(:,s)));
                            for hmm = 1:size(mapping.chromHMM,1) %%%%
                                if count_gene == 1,
%                                     count_features = count_features + 2;
                                    count_features = count_features + 1;
                                    features.names = [features.names; {[score_labels.operations{3} '_HMM:' ...
                                        mapping.chromHMM{hmm} '_' score_labels.regions{r}]}];
%                                         {[score_labels.operations{4} '_HMM:' ...
%                                         mapping.chromHMM{hmm} '_' score_labels.regions{r}]}];
                                end
                                
                                % # of rare variants for corresponding chromatin states
                                temp_idx.state = find(target_data.features(:,s) == hmm);
                                if ln == 0, % No chromatin states
%                                     temp_features = [temp_features NaN NaN];
                                    temp_features = [temp_features NaN];
                                elseif (ln > 0) && (length(temp_idx.state) > 0),
                                    temp_features = [temp_features sum(target_data.nvar(temp_idx.state))]; ...
%                                         length(temp_idx.state)/ln];
                                elseif (ln > 0) && (length(temp_idx.state) == 0),
%                                     temp_features = [temp_features 0 0];
                                    temp_features = [temp_features 0];
                                end
                            end %%%%
                        case 'DistTSS',
                            if count_gene == 1,
                                count_features = count_features + 1;
                                features.names = [features.names; {[score_labels.operations{2} '_' ...
                                    target_names{s} '_' score_labels.regions{r}]}];
                                features.region(r).pos_fnames = [features.region(r).pos_fnames ...
                                    {[score_labels.operations{2} '_' target_names{s} '_' score_labels.regions{r}]}];
                            end
                            [y,indices] = nanmin(target_data.features(:,s));
                            temp_features = [temp_features y];
                            temp_pos = [temp_pos target_data.pos(indices,:)];
                        otherwise,
                            if count_gene == 1,
                                count_features = count_features + 1;
                                features.names = [features.names; {[score_labels.operations{1} '_' ...
                                    target_names{s} '_' score_labels.regions{r}]}];
                                features.region(r).pos_fnames = [features.region(r).pos_fnames ...
                                    {[score_labels.operations{1} '_' target_names{s} '_' score_labels.regions{r}]}];
                            end
                            [y,indices] = nanmax(target_data.features(:,s));
                            temp_features = [temp_features y];
                            temp_pos = [temp_pos target_data.pos(indices,:)];
                    end
                end %%% for s = 1:length(target_names) %%%
                
                if count_gene == 1,
                    count_features = count_features + 1;
                    features.names = [features.names; {[score_labels.operations{3} '_all_' score_labels.regions{r}]}];
                end
                
                temp_features = [temp_features sum(target_data.nvar)]; % # of rare variants
                features.region(r).values = [features.region(r).values; temp_features];
                features.region(r).pos_oper = [features.region(r).pos_oper; temp_pos];
                
                if floor(g/10000) ~= factor,
                    factor = factor + 1;
                    disp([' *** ' score_labels.regions{r} ', g = ' num2str(g) ', # rv = ' num2str(count_gene) ' *** ']);
                end
            end
        end %% for g = 1:length(anno.idx_all) %%
        
        %% remove redundant phylop from Xin
        idx1_phylop = strmatch('max_phylop_10kb',features.region.pos_fnames,'exact');
        features.region(r).pos_fnames = features.region(r).pos_fnames';
        features.region(r).pos_fnames(idx1_phylop) = [];
        features.region(r).pos_fnames = features.region(r).pos_fnames';
        
        features.region(r).pos_oper = features.region(r).pos_oper';
        features.region(r).pos_oper(idx1_phylop,:) = [];
        features.region(r).pos_oper = features.region(r).pos_oper';
        
        idx2_phylop = strmatch('max_phylop_10kb',features.names,'exact');
        features.region(r).values = features.region(r).values';
        features.region(r).values(idx2_phylop,:) = [];
        features.region(r).values = features.region(r).values';
        
        features.names(idx2_phylop) = [];
    end % for r = 1:length(score_labels.regions) %
    time_stamp = toc/60
    disp([' ** ' num2str(ii) 'th indiv was done ... ** ']);    
    
    save(sprintf('%s/RIVER/data/score/indiv/%s.%d.score.gene.mat',dirNames.RAREVARDIR,id_ind,ii),'features');    
end
disp([' === ' num2str(time_stamp) ', # of features = ' num2str(length(features.names)) ' === ']);


