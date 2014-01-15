function [params] = convert_bam(n, params)
%CONVERT_BAM convert bam/sam file to sam/bam file
%   if sam file is passed in, converts to bam
%   if bam file is passed in, converts to sam
%   inputs:
%       file_in - sam/bam file to convert
%       file_out - sam/bam file to convert to

%% Accepted_hits
% extract name of file
[path, name, ext] = fileparts(char(params.accepted_bam_files{n}));
f_out = [path, '/', name, '.sam'];


if ~exist(f_out, 'file'),
    % if sam in and bam out
    if strcmp(params.accepted_bam_files{n}(end-3:end), '.bam') && strcmp(f_out(end-3:end), '.sam')
        convert_call = ['samtools view -h -o ', f_out, ' ', params.accepted_bam_files{n}];
        
        % if bam in and sam out
    elseif strcmp(params.accepted_bam_files{n}(end-3:end), '.sam') && strcmp(f_out(end-3:end), '.bam')
        convert_call = ['samtools view -hbS -o ', f_out, ' ', params.accepted_bam_files{n}];
    end
    
    unix(convert_call);
    
    % save .sam files in parameters
    params.accepted_sam_files{n} = f_out;
    
else
    [path, name, ext] = fileparts(params.align{n});
    samplename = strrep(path, params.outputdir, '');
    disp(['convert bam: ', samplename, ' accepted_sam file already exist']);
    % save .sam files in parameters
    params.accepted_sam_files{n} = f_out;
end
%% Unmapped_hits
% extract name of file
[upath, uname, uext] = fileparts(char(params.unmapped_bam_files{n}));
uf_out = [upath, '/', uname, '.sam'];

if ~exist(uf_out),
    % if sam in and bam out
    if strcmp(params.unmapped_bam_files{n}(end-3:end), '.bam') && strcmp(uf_out(end-3:end), '.sam')
        convert_call = ['samtools view -h -o ', uf_out, ' ', params.unmapped_bam_files{n}];
        
        % if bam in and sam out
    elseif strcmp(params.unmapped_bam_files{n}(end-3:end), '.sam') && strcmp(uf_out(end-3:end), '.bam')
        convert_call = ['samtools view -hbS -o ', uf_out, ' ', params.unmapped_bam_files{n}];
    end
    
    unix(convert_call);
    
    % save .sam files in parameters
    params.unmapped_sam_files{n} = uf_out;
    
else
    disp(['convert bam: ', samplename, ' unmapped_sam file already exist']); 
    % save .sam files in parameters
    params.unmapped_sam_files{n} = uf_out;
end
end

