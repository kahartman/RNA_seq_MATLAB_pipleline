function [chrom_names, gene_names, chrom, starts, stops, strand] = parse_gtf(params)
%PARSE_GTF extract relevant information from GTF File
%   Detailed explanation goes here

disp('Importing GTF file')

f_in = fopen(params.gtf_file);

tic;

format = '%s %s %s %d %d %s %s %s %s';

imported_txt = textscan(f_in, format, 'delimiter', '\t');

chrom_names = unique(imported_txt{1});

for i=1:numel(imported_txt{9})
    attributes = textscan(imported_txt{9}{i}, '%s%s%s%s%s%s%s%s%s%s');
    all_names{i} = attributes{2}{1};
end


disp('parsing GTF for annotations')

[gene_names, inds, groupings] = unique(all_names);
all_gene_inds = unique(groupings);

for i = 1:numel(all_gene_inds)
    gene_inds = groupings==all_gene_inds(i);
    starts(i) = min(imported_txt{4}(gene_inds));
    stops(i) = max(imported_txt{5}(gene_inds));
    
    chrom(i) = imported_txt{1}(inds(i));
    strand(i) = imported_txt{7}(inds(i));
end

save(params.mat_file, ...
     'chrom_names', 'gene_names', 'chrom', 'starts', 'stops', 'strand');

toc;


end

