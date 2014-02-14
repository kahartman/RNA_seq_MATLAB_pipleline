function [params] = runtophat(n, params)
% Run Tophat

% create tophat call for paired end or single end reads
if params.paired_ends
    base_name = params.fastq_left_files{n}(1:end-8);
    outfile = [params.sam_dir, base_name];
    
    run_tophat = ['nice ', params.tophat_dir, 'tophat ', ...
        strjoin(params.tophat_options), ' ', ...
        '--transcriptome-index ', params.gtf_file, ' ', ...
        '-o ', outfile, ' ', ...
        params.genome, ' ', ...
        params.fastq_dir, params.fastq_left_files{n} , ' ', ...
        params.fastq_dir, params.fastq_right_files{n} ...
        ];
else
    base_name = params.fastq_files{n}(1:end-6);
    outfile = [params.sam_dir, base_name];
    
    run_tophat = ['nice ', params.tophat_dir, 'tophat ', ...
        strjoin(params.tophat_options), ' ', ...
        '--transcriptome-index ', params.gtf_file, ' ', ...
        '-o ', outfile, ' ', ...
        params.genome, ' ', ...
        params.fastq_dir, params.fastq_files{n} ...
        ];
end



if ~exist(outfile, 'file'),
    % execute tophat command line call if file isn't already aligned
    system(run_tophat)
else
    disp(['tophat: ', outfile, ' already exists']);
end

% create cell of .bam files
params.bam_files{n} = [outfile, '/accepted_hits.bam'];


% convert .bam to .sam file
for n = 1:numel(params.bam_files),
    convert_bam(params.bam_files{n}, [params.sam_dir, base_name, '.sam']);
end;

end

