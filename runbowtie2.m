function [params] = runbowtie2(n, params )
%CALL_BOWTIE2 allign fastq_file using bowtie2 with specified parameters

if params.paired_ends
    % if reads are paired ends collect base names and create bowtie2 call
    sam_file = [params.fastq_left_files{n}(1:end-8), '.sam'];
    
    bowtie2_call = ['nice ', params.bowtie2_dir, 'bowtie2', ...
        ' ', strjoin(params.bowtie2_options, ' '), ' ', ...
        ' -x ', params.genome, ...
        ' -1 ', [params.fastq_dir, params.fastq_left_files{n}], ...
        ' -2 ', [params.fastq_dir, params.fastq_right_files{n}], ...
        ' -S ', [params.sam_dir, sam_file] ...
        ];
    
else
    % bowtie2 call for single end reads
    sam_file = [params.fastq_files{n}(1:end-6), '.sam'];
    
    bowtie2_call = ['nice ', params.bowtie2_dir, 'bowtie2', ...
        ' ', strjoin(params.bowtie2_options, ' '), ' ', ...
        ' -x ', params.genome, ...
        ' -U ', [params.fastq_dir, params.fastq_files{n}], ...
        ' -S ', [params.sam_dir, sam_file] ...
        ];
end

if ~exist([params.sam_dir sam_file], 'file')
    
    % execute bowtie2 call
    disp('Running Bowtie2')
    tic;
    unix(bowtie2_call);
    time = toc;
    disp([sam_file, ' complete, ', 'time:', num2str(time)])
    
else
    % skip alignment if already aligned
    disp(['bowtie2: ', sam_file, ' already exist']);
end

% create cell of .bam files
params.sam_files{n} = sam_file;



end
