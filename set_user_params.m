function params = set_user_params()
%SET_PIPELINE_PARAMS user inputed parameters
%   These parameters may need to be changed from run to run or user to user


% input parameters
%params.fastq_dir = '/home/david/RNAseq_files/MiSeq_Mar21_2014/fastqs/';
%params.fastq_dir = '/home/david/RNAseq_files/tissue_test/fastqs/';
params.fastq_dir = '/home/david/RNAseq_files/Brn2Data/fastqs/';
%params.fastq_dir = '/home/david/RNAseq_files/tissues/fastqs/';

% collect single or paired end reads
params.paired_ends = 0;


% sam output path
%params.sam_dir = '/home/david/RNAseq_files/MiSeq_Mar21_2014/alignments/';
%params.sam_dir = '/home/david/RNAseq_files/tissue_test/alignments/';
params.sam_dir = '/home/david/RNAseq_files/Brn2Data/alignment/';
%params.sam_dir = '/home/david/RNAseq_files/tissues/alignments/';

% choose to run tophat (1) or bowtie (0)
params.tophat = 0;
% align to transcriptome(1) or genome (0)
params.transcriptome = 1;
% number of processors to use for alignments
params.processors = 8;


params.genome = '/home/david/RNAseq_files/reference_genomes/UCSC/mm9/Transcriptome/transcriptome';
params.gtf_file = '/home/david/RNAseq_files/reference_genomes/UCSC/mm9/Annotation/Genes/genes.gtf';



if params.tophat == 0
    % bowtie2 parameters
    params.bowtie2_dir = '';
    params.bowtie2_options = { '-D', '25', '-R', '3', '-N', '1', '-L', '20',...
        '-i', 'S,1,0.50', '--local',...
        '-p', num2str(params.processors) ...
        };
else
    % tophat parameters
    params.tophat_dir = '';    
    params.tophat_options = { '-p', num2str(params.processors), ...
        %'--no-coverage-search', '-N', '4', ...
        %'--read-edit-dist', '4', '--b2-N', '1' ...
        };


%%%%%%%%%%%%%%%    params.transcriptome =  '/home/graham/Projects/reference_genomes/UCSC/mm9/Annotation/Genes/genes';

end

% filename for completed params data struct (ending should be '.mat')
%params.output_filename = '/home/david/RNAseq_files/MiSeq_Mar21_2014/Mar21_2014_data_allgenes.mat';
%params.output_filename = '/home/david/RNAseq_files/tissue_test/tissue_test.mat';
params.output_filename = '/home/david/RNAseq_files/Brn2Data/Brn2_data_all.mat';
%params.output_filename = '/home/david/RNAseq_files/tissues/tissue_data.mat';


end

