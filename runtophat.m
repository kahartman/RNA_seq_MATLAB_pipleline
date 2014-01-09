function [params] = runtophat(n, params)
% params = runtophat(fastq_files_Q, params)

%% Run Tophat
% modify output file name
[path, name, ext] = fileparts(char(params.align{n}));
% find sample folder name
qualscore = char(strcat('_Q', num2str(params.qual)));
split = strsplit(name, qualscore);
fastq_name = split{1};
outfile = [params.outputdir, fastq_name, '/', name, '/accepted_hits.bam'];

if ~exist(outfile, 'file'),
    % run tophat
    run_tophat = [params.tophat_dir, 'tophat ', ...
        strjoin(params.tophat_options), ' ', ...
        '--transcriptome-index ', params.gtf_file, ' ', ...
        '-o ', params.sam_dir, fastq_name, '/', name, '/ ', ...
        params.genome, ' ', ...
        params.outputdir, fastq_name, '/', name, ext, ...
        ];
    
    system(run_tophat)
    
    % create cell of .bam files
    params.accepted_bam_files{n} = [params.outputdir, fastq_name, '/', name, '/accepted_hits.bam'];
    params.unmapped_bam_files{n} = [params.outputdir, fastq_name, '/', name, '/unmapped.bam'];
    
else
    disp(['tophat: ', fastq_name, ' already exist']);
    % create cell of .bam files
    params.accepted_bam_files{n} = [params.outputdir, fastq_name, '/', name, '/accepted_hits.bam'];
    params.unmapped_bam_files{n} = [params.outputdir, fastq_name, '/', name, '/unmapped.bam'];
    
end

%% convert .bam to .sam file
for n = 1:numel(params.accepted_bam_files),
   [params] = convert_bam(n, params);
end;
