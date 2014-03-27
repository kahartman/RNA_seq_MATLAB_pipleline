function params = set_pipeline_params()


params = set_user_params();


split_name = strsplit(params.gtf_file, '.');
params.mat_file = [split_name{1}, '.mat'];
if exist(params.mat_file, 'file')
    disp('Loading pre-existing .mat annotation file')
    load(params.mat_file);
    params.chrom_names = chrom_names;
    params.gene_names = gene_names;
    params.chroms = chrom;
    params.strands = strand;
else
    [params.chrom_names, params.gene_names, params.chroms, starts, stops, params.strands] ...
        = parse_gtf(params);
end

params.lengths = stops - starts;

params.chrom_num = numel(params.chrom_names);

params.gene_annotation_hash = containers.Map(params.gene_names, ...
                                            1:numel(params.gene_names));


% setup gene start/stop data structures
params.chrom_hash = containers.Map(params.chrom_names, ...
                                   1:numel(params.chrom_names));
params.rev_chrom_hash = containers.Map(1:numel(params.chrom_names), ...
                                       params.chrom_names);
                               
[~,~,c] = unique(params.chroms);                               
max_genes_on_chrom = sum(c == mode(c));

params.starts = zeros(params.chrom_num, max_genes_on_chrom);
params.stops = zeros(params.chrom_num, max_genes_on_chrom);

for i=1:params.chrom_num
    curr_chrom_inds = strcmp(params.rev_chrom_hash(i), params.chroms);

    params.starts(i,1:numel(starts(curr_chrom_inds))) = starts(curr_chrom_inds);
    params.stops(i,1:numel(stops(curr_chrom_inds))) = stops(curr_chrom_inds);
    
end


params.allsamples = [];




end

