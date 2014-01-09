%% import user defined parameters

close all; clear all; clc;
params = set_pipeline_params();

%% grab .fastq files
%params = grabfiles(params);

%% remove low quality reads
%for n = 1:numel(params.fastq_files), 
%    params = qualityfilter(n, params);
%end
%% run fastQC
%params = fastqc(params);

%% trim low quality bp from front/back of reads
%for n = 1:numel(params.align),
%    params = trimseq(n, params);
%end

%% align reads with tophat or bowtie
for n = 1:numel(params.fastq_files)
    if params.tophat
        params = runtophat(n, params);
    else
        params = runbowtie2(n, params);
    end
end

%% read accepted_hits.sam files
for n = 1:numel(params.sam_files),
   params = genecounts(n, params);
end

