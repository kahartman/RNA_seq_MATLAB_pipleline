function [aligned_count, counts] = hash_sam_gene_align(read_starts, read_chroms, params)

% read sam file
[start, chrom] = read_sam([params.sam_dir, params.sam_files{n}]);



% Align to genes

counts = zeros(size(params.starts));

aligned_count = 0;
disp('mapping reads to genes');

% loop over all chroms
for curr_chrom = 1:params.chrom_num
    
    % identify reads mapping to the current chrom
    chrom_inds = strcmp(params.chrom_names{curr_chrom}, read_chroms);
    
    % loop over all genes on a chrom
    for curr_gene = 1:numel(params.starts(curr_chrom)>0),
        
        % find reads mapping to current gene (on current chrom)
        hits = ((read_starts(chrom_inds) > params.starts(curr_chrom, curr_gene)) & ...
                (read_starts(chrom_inds) < params.stops(curr_chrom, curr_gene)));
        
        % if the read aligns to a gene, count it within the gout matrix
        if sum(hits)>0
            counts(curr_chrom, curr_gene) = counts(curr_chrom, curr_gene) + sum(hits);
            aligned_count = aligned_count + sum(hits);
        end;
        
%         if mod(read_number, 100) == 0,
%             disp([num2str(read_number), ' reads mapped']);
%         end;
    end;
end;

params.allsamples = [params.allsamples, counts(params.starts>0)];
end

