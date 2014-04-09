%% import user defined parameters

close all force; clear all; clc;
cd('/home/david/MATLAB/RNAseq_pipeline/');
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

    matlabpool(3);
    counts_hash = cell(1,params.file_num);

    time = zeros(1, params.file_num);
    parfor n = 1:params.file_num
        time(n) = tic;
        counts_hash{n} = sam_gene_align(n, params);
        disp([num2str(n), ':', num2str(numel(params.sam_files))])
        toc(time(n))
    end

    params.counts_hash = counts_hash;
    
    
    all_counts = zeros(params.gene_num, numel(params.counts_hash));
    
    for gene_ind = 1:params.gene_num    
        for exp_ind = 1:numel(params.counts_hash)
                all_counts(gene_ind, exp_ind) = params.counts_hash{exp_ind}(params.gene_names{gene_ind});
        end
    end
    
    params.allsamples = all_counts;
    save(params.output_filename, 'params');
    matlabpool close

else 
    load(params.output_filename)
end

disp('params loaded in')



