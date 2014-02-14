function convert_bam( f_in, f_out )
%CONVERT_BAM convert bam/sam file to sam/bam file
%   if sam file is passed in, converts to bam
%   if bam file is passed in, converts to sam
%   inputs:
%       file_in - sam/bam file to convert
%       file_out - sam/bam file to convert to

% if sam in and bam out
if strcmp(f_in(end-3:end), '.bam') && strcmp(f_out(end-3:end), '.sam')
    convert_call = ['samtools view -h -o ', f_out, ' ', f_in];

% if bam in and sam out
elseif strcmp(f_in(end-3:end), '.sam') && strcmp(f_out(end-3:end), '.bam')
    convert_call = ['samtools view -hbS -o ', f_out, ' ', f_in];  
end

% execute samtools call at command line
unix(convert_call);

end

