function params = genecounts(n, params)


% read sam file
[start, chrom] = read_sam([params.sam_dir, params.sam_files{n}]);
% run gene alignment
[aligned_counts, gout_vec] = hash_sam_gene_align(start, chrom, params);

% combine sample to others in project
params.allsamples = [params.allsamples, gout_vec];

