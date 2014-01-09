function [params] = runbowtie2(n, params )
%CALL_BOWTIE2 allign fastq_file using bowtie2 with specified parameters

sam_file = [params.fastq_files{n}(1:end-6), '.sam'];
tic;
if ~exist([params.sam_dir sam_file], 'file')
    disp('Running Bowtie2')
    bowtie2_call = [params.bowtie2_dir, 'bowtie2', ...
        ' ', strjoin(params.bowtie2_options, ' '), ' ', ...
        ' -x ', params.genome, ...
        ' -U ', [params.fastq_dir, params.fastq_files{n}], ...
        ' -S ', [params.sam_dir, sam_file] ...
        ];
    disp(bowtie2_call)
    
    
    unix(bowtie2_call);
    
    % create cell of .bam files
    params.sam_files{n} = sam_file;
    
else
    disp(['bowtie2: ', sam_file, ' already exist']);
    % create cell of .bam files
    params.sam_files{n} = sam_file;
    
end
toc;
end
