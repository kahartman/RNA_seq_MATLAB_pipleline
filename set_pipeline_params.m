function params = set_pipeline_params()
%SET_PIPELINE_PARAMS user inputed parameters
%   These parameters may need to be changed from run to run or user to user


% input parameters
params.fastq_dir = '/home/graham/MATLAB/RNAseq_pipeline/fastqs/';

fastqs = dir([params.fastq_dir, '*.fastq']);
params.fastq_files = extractfield(fastqs, 'name');



% output parameters
params.sam_dir = '/home/graham/MATLAB/RNAseq_pipeline/alignments/';

% % fastx parameters
% params.fastxdir = '/home/thomsonlab/Bowtie/fastx_toolkit-0.0.12/';
% params.qual = '30';             % Minimum quality score to keep- fastq_quality_filter
% params.pct = '75';              % Minimum percent of bases that must have -q quality
% params.sng = '33';              % 33 for Sanger sequences
% params.gatc_content_threshold = 33;
% params.mean_quality_threshold = 30;

% % fastqc parameters
% params.fastqc_txt_files = {};
% params.fastqcdir = '/home/thomsonlab/Bowtie/FastQC/';

% choose to run tophat (1) or bowtie (0)
params.tophat = 0;
params.processors = 4;
params.genome = '/home/graham/Projects/reference_genomes/Gloria/gloria_sk1/gloria_sk1';
%params.genome = '/home/graham/Projects/reference_genomes/UCSC/mm9/Sequence/Bowtie2Index/genome';
params.gtf_file = '/home/graham/Projects/reference_genomes/Gloria/genes/genes_only.gff';
%params.gtf_file = '/home/graham/Projects/reference_genomes/UCSC/mm9/Annotation/Genes/genes.gtf';
if params.tophat == 0
    % bowtie2 parameters
    params.bowtie2_dir = '/home/graham/Software/bowtie2-2.1.0/';
    params.bowtie2_options = { '-D', '25', '-R', '3', '-N', '1', '-L', '20',...
        '-i', 'S,1,0.50', '--local',...
        '-p', num2str(params.processors) ...
        };
else
    % tophat parameters
    params.tophat_dir = '/home/graham/Bowtie/Tophat/';
    params.transcriptome =  '/home/graham/Projects/reference_genomes/UCSC/mm9/Annotation/Genes/genes';
    params.tophat_options = { '-p', num2str(params.processors), ...
        '--no-coverage-search', '-N', '4', ...
        '--read-edit-dist', '4', '--b2-N', '1' ...
        };
end



% annotation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

params.chrom_num = numel(params.chrom_names);


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
    [sorted_starts, sort_inds] = sort(starts(curr_chrom_inds));

    sorted_stops = stops(curr_chrom_inds);
    sorted_stops = sorted_stops(sort_inds);
    
    params.starts(i,1:numel(sorted_starts)) = sorted_starts;
    params.stops(i,1:numel(sorted_stops)) = sorted_stops;
    
end


params.allsamples = [];




end

% 
% 
% % Jade's Parameters
% 
% % pipeline timer
% tic1 = tic;
% time = [];
% 
% % input parameters
% params.fastq_dir = '/media/Thomson/raw_data/';
% 
% % output parameters
% params.outputdir =  '/media/Thomson/Josh_index1_newalign/';
% 
% % fastx parameters
% params.fastxdir = '/home/thomsonlab/Bowtie/fastx_toolkit-0.0.12/';
% params.qual = '30';             % Minimum quality score to keep- fastq_quality_filter
% params.pct = '75';              % Minimum percent of bases that must have -q quality
% params.sng = '33';              % 33 for Sanger sequences
% params.gatc_content_threshold = 33;
% params.mean_quality_threshold = 30;
% 
% % fastqc parameters
% params.fastqc_txt_files = {};
% params.fastqcdir = '/home/thomsonlab/Bowtie/FastQC/';
% 
% % bowtie2 parameters
% params.processors = 4;
% params.genome = '/home/thomsonlab/Bowtie/bowtie2-2.1.0/mm9/mm9';
% params.bowtie2_dir = '/home/thomsonlab/Bowtie/bowtie2-2.1.0/';
% params.bowtie2_options = { '-D', '25', '-R', '3', '-N', '1', '-L', '20',...
%                            '-i', 'S,1,0.50', '--local',...
%                            '-p', num2str(params.processors) ...
%                           };
% 
% 
% % % tophat parameters
% % params.tophat_dir = '/home/thomsonlab/Bowtie/Tophat/';
% % params.genome = '/home/thomsonlab/Bowtie/bowtie2-2.1.0/mm9/mm9';
% % params.transcriptome = '~/Bowtie/bowtie2-2.1.0/mm9/genes';
% % params.tophat_processors = 8;
% % params.tophat_options = { '-p', num2str(params.tophat_processors), ...
% %                         '--no-coverage-search', '-N', '4', ...
% %                         '--read-edit-dist', '4', '--b2-N', '1' ...
% %                         };
% 
% % alignment parameters
% params.project_name = 'Index1file';
% params.genefile = '/home/thomsonlab/Bowtie/bowtie2-2.1.0/mm9/mm9cord.csv';
% params.allsamples = [];
% 
