%% import user defined parameters

close all force; clear all; clc;
cd('~/MATLAB/RNAseq_pipeline/rna_seq_pipe/');
params = set_pipeline_params();


%% align reads with tophat or bowtie
for n = 1:params.file_num
    if params.tophat
        params = runtophat(n, params);
    else
        params = runbowtie2(n, params);
    end
end

%% read accepted_hits.sam files
if ~exist(params.output_filename)
    for n = 1:numel(params.sam_files)
        %for n = length(params.counts_hash):numel(params.sam_files)
        tic;
        params = sam_gene_align(n, params);
        disp([num2str(n), ':', num2str(numel(params.sam_files))])
        toc;
    end
    
    gene_names = [];
    for i=1:numel(params.counts_hash)
        gene_names = [gene_names, params.counts_hash{i}.keys()];
    end
    
    gene_names = unique(gene_names);
    master_gene_hash = containers.Map(gene_names, 1:numel(gene_names));
    
    all_counts = zeros(numel(gene_names), numel(params.counts_hash));
    
    for gene_ind = 1:numel(gene_names)
        
        %    if params.gene_annotation_hash.isKey(gene_names{gene_ind})
        %        lengths(gene_ind) = params.lengths(params.gene_annotation_hash(gene_names{gene_ind}));
        
        for exp_ind = 1:numel(params.counts_hash)
            if params.counts_hash{exp_ind}.isKey(gene_names(gene_ind))
                all_counts(gene_ind, exp_ind) = params.counts_hash{exp_ind}(gene_names{gene_ind});
            else
                all_counts(gene_ind, exp_ind) = 0;
            end
        end
        
        %    end
    end
    
    params.allsamples = all_counts;
    params.gene_names_ordered = gene_names;
    save('params.output_filename', 'params');
else 
    load('params.output_filename')
    
end
disp('params loaded in')



