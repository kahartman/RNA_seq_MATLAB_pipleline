function [params] = sam_gene_align(n, params)



% read sam file
[read_starts, read_chroms] = read_sam([params.sam_dir, params.sam_files{n}]);

% remove unaligned reads
unaligned_inds = strcmp('*', read_chroms);
read_starts(unaligned_inds) = [];
read_chroms(unaligned_inds) = [];

% remove duplicate reads
[single_copy_starts, ~, repeat_starts] = unique(read_starts);
[single_copy_chroms,~, repeat_chroms] = unique(read_chroms);
replicate_reads = repeat_starts == repeat_chroms;
read_starts(replicate_reads) = [];
read_chroms(replicate_reads) = [];
read_starts = [read_starts, single_copy_starts];
read_chroms = [read_chroms; single_copy_chroms]';


if params.transcriptome
    
    % each experiment has a hash tabel from genes to number of reads
    % if gene is not in hash, it is becuase it has zero reads mapping in
    % the experiment
    
    [gene_names, ~, gene_inds] = unique(read_chroms);
    for i=1:max(gene_inds)
        gene_counts(i) = sum(gene_inds == i);
    end
    
    params.counts_hash{n} = containers.Map(gene_names, gene_counts);
    
else
    % Align to genes
    counts = zeros(size(params.starts));
    
    % loop over all chroms
    for curr_chrom = 1:params.chrom_num
        
        % identify reads mapping to the current chrom
        chrom_inds = strcmp(params.chrom_names{curr_chrom}, read_chroms);
        
        % loop over all genes on a chrom
        for curr_gene = 1:numel(params.starts(curr_chrom,:)>0),
            
            % find reads mapping to current gene (on current chrom)
            hits = ((read_starts(chrom_inds) >= params.starts(curr_chrom, curr_gene)) & ...
                (read_starts(chrom_inds) <= params.stops(curr_chrom, curr_gene)));
            
            % if the read aligns to a gene, count it within the gout matrix
            if sum(hits)>0
                counts(curr_chrom, curr_gene) = counts(curr_chrom, curr_gene) + sum(hits);
            end;
            
        end;
    end;
    params.allsamples = [params.allsamples, counts(params.starts>0)];
    
    disp(params.fastq_files{n})
    disp([num2str(numel(read_starts)), ' reads']);
    disp([num2str(sum(sum(counts))), ' mapped to coding regions']);
    disp('')
end

end

