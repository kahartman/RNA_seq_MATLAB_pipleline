function params = set_user_params()
%SET_PIPELINE_PARAMS user inputed parameters
%   These parameters may need to be changed from run to run or user to user


% input parameters
%params.fastq_dir = '/home/graham/MATLAB/RNAseq_pipeline/fastqs/';
%params.fastq_dir = '/media/newhd/Projects/amit_data/fastqs/';
params.fastq_dir = '/media/newhd/Projects/regev_data/fastqs/';
%params.fastq_dir = '/media/newhd/Projects/tissues/fastqs/';
%params.fastq_dir = '/media/newhd/Projects/brar/fastqs/';

% collect single or paired end reads
params.paired_ends = 1;
if params.paired_ends
    fastqs_left = dir([params.fastq_dir, '*.fastq_1']);
    fastqs_right = dir([params.fastq_dir, '*.fastq_2']);
    params.fastq_left_files = extractfield(fastqs_left, 'name');
    params.fastq_right_files = extractfield(fastqs_right, 'name');
    params.file_num = numel(fastqs_left);
else
    fastqs = dir([params.fastq_dir, '*.fastq']);
    params.fastq_files = extractfield(fastqs, 'name');
    params.file_num = numel(fastqs);
end



% output parameters
%params.sam_dir = '/home/graham/MATLAB/RNAseq_pipeline/alignments/';
%params.sam_dir = '/media/newhd/Projects/amit_data/alignments/';
params.sam_dir = '/media/newhd/Projects/regev_data/alignments/';
%params.sam_dir = '/media/newhd/Projects/tissues/alignments/';
%params.sam_dir = '/media/newhd/Projects/brar/alignments/';



% choose to run tophat (1) or bowtie (0)
params.tophat = 0;
% align to transcriptome(1) or genome (0)
params.transcriptome = 1;
% number of processors to use for alignments
params.processors = 4;

%params.genome = '/home/graham/Projects/reference_genomes/UCSC/mm9/Sequence/Bowtie2Index/genome';
params.genome = '/media/newhd/Projects/reference_genomes/UCSC/mm9/Transcriptome/transcriptome';
params.gtf_file = '/media/newhd/Projects/reference_genomes/UCSC/mm9/Annotation/Genes/genes.gtf';
%params.genome = '/media/newhd/Projects/reference_genomes/Gloria/genes';
%params.gtf_file = '/media/newhd/Projects/reference_genomes/Gloria/genes_only.gff';

if params.tophat == 0
    % bowtie2 parameters
    params.bowtie2_dir = '/home/graham/Software/bowtie2-2.1.0/';
    params.bowtie2_options = { '-D', '25', '-R', '3', '-N', '1', '-L', '20',...
        '-i', 'S,1,0.50', '--local',...
        '-p', num2str(params.processors) ...
        };
else
    % tophat parameters
    params.tophat_dir = '/home/graham/Software/tophat-2.0.9.Linux_x86_64/';    
    params.tophat_options = { '-p', num2str(params.processors), ...
        %'--no-coverage-search', '-N', '4', ...
        %'--read-edit-dist', '4', '--b2-N', '1' ...
        };

%    params.transcriptome =  '/home/graham/Projects/reference_genomes/UCSC/mm9/Annotation/Genes/genes';

end


end
