%%%%% Assign relevant directories into your own directories 
dirNames.anno = '${RAREVARDIR}/reference';
dirNames.region = '${RAREVARDIR}/preprocessing/rvsite';
%%%%%

region_tss.labels = '10k';
region_tss.dist = 1e3*10;
region_tss.fid = fopen(sprintf('%s/region.tss%s.txt',dirNames.region,region_tss.labels),'w');
load(sprintf('%s/annoGTEx_112415.mat',dirNames.anno));

factor = 0;
for g = 1:size(anno.gene_id,1)
    if anno.chr(g) < 23, % remove X,Y,MT
        switch anno.gene_id{g,1},
            case '+',
                idx_tss = 1;
            case '-',
                idx_tss = 2;
        end
                
        for i = 1:length(region_tss.labels)
            % chr | start | stop | gencodeidx | geneName | TSS
            fprintf(region_tss.fid(i),'%d\t%d\t%d\t%d\t%s\t%d\n',anno.chr(g),max(1,anno.gene_pos(g,idx_tss)-region_tss.dist(i)),anno.gene_pos(g,idx_tss)+region_tss.dist(i),anno.idx_all(g),anno.gene_id{g,2},anno.gene_pos(g,idx_tss));
        end
    else
        continue
    end
    if floor(g/1000) ~= factor,
        factor = factor + 1;
        disp([' === g: ' num2str(g) ' === ']);
    end
end
fclose(region_tss.fid);

    
                
                
                
    
