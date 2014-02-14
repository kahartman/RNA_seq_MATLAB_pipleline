%% import user defined parameters

close all; clear all; clc;
params = set_pipeline_params();

%% grab .fastq files
%params = grabfiles(params);

%% remove low quality reads
%for n = 1:numel(params.fastq_files), 
%    params = qualityfilter(n, params);
%end
%% run fastQC
%params = fastqc(params);

%% trim low quality bp from front/back of reads
%for n = 1:numel(params.align),
%    params = trimseq(n, params);
%end

%% align reads with tophat or bowtie
for n = 1:params.file_num
    if params.tophat
        params = runtophat(n, params);
    else
        params = runbowtie2(n, params);
    end
end

%% read accepted_hits.sam files
for n = 1:numel(params.sam_files)
%for n = length(params.counts_hash):numel(params.sam_files)
   tic;
   params = sam_gene_align(n, params);
   disp([num2str(n), ':', num2str(numel(params.sam_files))])
   toc;
end

%save('test_params.mat', 'params');

%%
gene_names = [];
for i=1:numel(params.counts_hash)
    gene_names = [gene_names, params.counts_hash{i}.keys()];
end

gene_names = unique(gene_names);
master_gene_hash = containers.Map(gene_names, 1:numel(gene_names));

all_counts = zeros(numel(gene_names), numel(params.counts_hash));

for gene_ind = 1:numel(gene_names)
    
%    if params.gene_annotation_hash.isKey(gene_names{gene_ind})
%        lengths(gene_ind) = params.lengths(params.gene_annotation_hash(gene_names{gene_ind}));
        
        for exp_ind = 1:numel(params.counts_hash)
            if params.counts_hash{exp_ind}.isKey(gene_names(gene_ind))
                all_counts(gene_ind, exp_ind) = params.counts_hash{exp_ind}(gene_names{gene_ind});
            else
                all_counts(gene_ind, exp_ind) = 0;
            end
        end
        
%    end
end

% %% calculating KL-divergence
% 
% figure;
% loglog(all_counts(:,1), all_counts(:,29), 'o', 'LineWidth' , .5)
% xlabel('Cell 1')
% ylabel('Cell 10')
% title('Similarity of gene expression between cells')
% 
% [a,b,c] = pca(all_counts'./repmat(sum(all_counts), size(all_counts,1),1)');
% figure;
% plot(cumsum(c)/sum(c), 'x-')
% xlabel('PCs')
% ylabel('% varience explained')
% label('tight')
% 
% [idx, centers] = kmeans(all_counts', 4, 'Distance', 'correlation');
% 
% col_1 = all_counts(:,1) + 1*ones(size(all_counts(:,1)));
% col_2 = all_counts(:,29) + 1*ones(size(all_counts(:,5)));
% 
% n = max(col_1, col_2);
% eps = min(col_1, col_2) ./ max(col_1, col_2);
% 
% % figure;
% % semilogx(n,eps, '.')
% % axis('tight')
% % xlabel('number of reads')
% % ylabel('epsilon')
% 
% kl = .5*(log(eps) + (1 + n)./eps + n.*(eps - 2) - 1);
% 
% % figure;
% % loglog(n, kl, '.')
% % axis('tight')
% % xlabel('number of reads')
% % ylabel('KL divergence')
% 
% figure;
% scatter(log10(n), eps, 10, log10(kl));
% line([log10(4054), log10(4054)], [0,1])
% axis('tight')
% xlabel('number of reads')
% ylabel('epsilon')
% title('Information from each gene')
% 
% figure;
% hist(eps, 20)
% xlabel('epsilon')
% ylabel('frequency')
% title('Distribution of epsilons')
% 
% 
% 
% 
% 
% %% Downsample data
% 
% sample_num_vec = 1e4:.5e5:2e6;
% 
% for sample_rate_ind = 1:numel(sample_num_vec) 
%     sample_num = sample_num_vec(sample_rate_ind);
%     for n = 1:numel(params.sam_files)
%         downsample{sample_rate_ind}(:,n) = sample_distribution(params.allsamples(:,n), sample_num);
%     end
% 
%     [W{sample_rate_ind}, H{sample_rate_ind}] = nnmf(downsample{sample_rate_ind}, 10, 'algorithm', 'mult');
% 
% end
% 
% 
% %%
% for i=1:numel(W)
%     W_aligned{i} = align_matrices(W{i}, W{end});
% end
% 
% %%
% 
% for i=1:7
%     mapobj = HeatMap(W_aligned{i}, 'symmetric', 0, 'standardize', 1, 'DisplayRange', 1);
%     Mh(i) = getframe(view(mapobj));
% end
% movie2avi(Mh, '/home/graham/Desktop/test.avi', 'fps', 2)








