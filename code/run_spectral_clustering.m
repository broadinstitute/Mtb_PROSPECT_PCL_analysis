% # Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, false);
%imatlab_export_fig('print-png')

% loocv inputs

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'number' % 'number' or 'list'

demo_loocv_number_cmpds = 1

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359','BRD-K87202646','BRD-K59853741', 'BRD-K27302037'} % Ciprofloxacin, Rifampin, Isoniazid, Q203, Thioacetazone

results_subdir_prefix = 'loocv_pcls/leave_out_cmpd_'

loocv_save_fig = true % .png files of spectral clustering input and output

loocv_save_out = true % save tabular and gmt files for each MOA separately with the output of spectral clsutering

unique_kabx_cmpds_tbl_path = '../results/kabx_pert_ids_tbl_for_loocv.txt'


% general inputs

moa_figdir = 'verify_moas/figures'

moa_gctdir = 'verify_moas/gcts'

corr_for_pcls_savepath = 'sGR_for_pcls_pearson_corr'

corr_rank_for_pcls_savepath = 'sGR_for_pcls_pearson_corr_rank'

moa_gmt_savepath = 'moas.gmt'

save_fig = true % .png files of spectral clustering input and output

save_out = true % save tabular and gmt files for each MOA separately with the output of spectral clsutering


% Spectral Clustering inputs

rng_seed = 0; % specified seed for initializing the random number generator, Matlab factory default is the Mersenne Twister generator with seed 0 (see spectral_clustering_for_pcls.m for additional information)

num_threads_for_multithreading = 4 % number of computational threads to use for multi-threading enabled Matlab functions including eigs; data was originally processed with 4 CPU cores or threads, this is set for reproducibility across CPU hardware

thrsh_rank = 20 % threshold for average pairwise rank of correlation across KABX to connect treatments as mutual nearest-neighbors

thrsh_factor = 1; % factor to multiply by thrsh_rank if dynamic_thrsh_per_moa is false; default is 1

dynamic_thrsh_per_moa = false % if true then threshold is round(log(size of MOA) * thrsh_rank), otherwise identical threshold for every MOA

k_type = 'k_med_gap_den' % eigengap heuristic to take for estimating number of K clusters: k_num_zero, k_num_zero_plus_one, k_med_gap_den, k_gap_den (see create_laplacian_matrix.m for additional information)

show_hclust = true % if true then apply hierarchical clustering to correlation matrix prior to running spectral clustering, whether or not this is performed can minimally impact final cluster results especially for larger MOAs with greater number of K clusters due to stochasticity (randomized initialization) of k-means++ clustering; original results did apply hierarchical clustering but this step is not required for the method to be successful

% outputs

if dynamic_thrsh_per_moa
    outdir_name = sprintf('clusters_spectral_clustering_thrsh_rank_le%dxlogsize_%s', thrsh_rank, k_type) % log(MOA size), i.e. the number of treatments/dsCGI profiles in the MOA
else
    outdir_name = sprintf('clusters_spectral_clustering_thrsh_rank_le%d_%s', thrsh_rank, k_type)
end

outdir = fullfile(wkdir, outdir_name)

moa_k_values_filename = 'k_values.txt'

clusters_tbl_incl_singletons_filename = 'clusters_spectral_clust_including_singletons.txt'

clusters_tbl_filename = 'clusters_spectral_clust.txt'

clusters_gmt_filename = 'clusters_spectral_clust.gmt'

disp(sprintf('Number of computational threads: %d', maxNumCompThreads));
maxNumCompThreads(num_threads_for_multithreading); % to reproduce results set num_threads_for_multithreading to 4 as it was originally processed with 4 available CPU cores 
disp(sprintf('Number of computational threads: %d', maxNumCompThreads));

switch k_type
case 'k_gap_den'
case 'k_med_gap_den'
case 'k_num_zero'
case 'k_num_zero_plus_one'
otherwise
    error('Unknown k_type; Only "k_gap_den", "k_med_gap_den", "k_num_zero_plus_one", or "k_num_zero" can be used')
end

% Ensure parallel pool is cleaned up at the end of the script or on error
cleanupObj = onCleanup(@() delete(gcp('nocreate')));

% Commented out/disabled parallel processing on Code Ocean and for setting random seed
% p = parpool(feature('numcores'),'IdleTimeout',120);

moas = parse_gmt(fullfile(wkdir,moa_gmt_savepath));

mk_cd_dir(outdir, false);
mk_cd_dir(wkdir, false);

% # Run Spectral Clustering (and Plot Eigenvalues) for Each MOA

spec = {};
out_lap = {};

out = {};
out_nonsingle = {};
gmt = {};

num_moas_represented = 0;


for ii = 1:numel(moas)
    loop_progress(ii,numel(moas),1)
    
    % Set MoA class
    moa_class = moas(ii).head;
    out_lap(ii).moa_class = moa_class;
    disp(moa_class)
    
    %[strrep(strrep(moa_class,' ','_'),'/','-'),'_corr_n*.gctx']
    
    % Load gcts
    g = glob(fullfile(wkdir, moa_gctdir,[strrep(strrep(moa_class,' ','_'),'/','-'),'_corr_n*.gctx']));
    
    try
        c = parse_gctx(g{1});
    catch ME
        disp(sprintf('Error: %s correlation gctx file not found', moa_class))
        disp(ME)
    end
    
    g = glob(fullfile(wkdir, moa_gctdir,[strrep(strrep(moa_class,' ','_'),'/','-'),'_corr_rank_n*.gctx']));
    
    try
        cr = parse_gctx(g{1});
    catch ME
        disp(sprintf('Error: %s correlation rank gctx file not found', moa_class))
        disp(ME)
    end
    
    moa_size = numel(cr.rid);
    log_moa_size = log(moa_size);
    
    % Set threshold
    if dynamic_thrsh_per_moa
        thrsh = round(log_moa_size * thrsh_rank);
    else
        thrsh = thrsh_rank*thrsh_factor;
    end
    
    out_lap(ii).moa_size = moa_size;
    out_lap(ii).log_moa_size = log_moa_size;
    out_lap(ii).thrsh_rank = thrsh_rank;
    out_lap(ii).thrsh_factor = thrsh_factor;
    out_lap(ii).thrsh = thrsh;
    
    % Run spectral clustering
    try
        [spec(ii).tmp_out_nonsingle,spec(ii).tmp_out,spec(ii).tmp_gmt,L,out_lap(ii).k,en,den,out_lap(ii).k_gap_den,out_lap(ii).k_med_gap_den,out_lap(ii).k_num_zero_plus_one,out_lap(ii).k_num_zero] = spectral_clustering_for_pcls(c,cr,moa_class,thrsh,k_type,outdir,save_out,save_fig,rng_seed,show_hclust);
        close all
        
        if size(spec(ii).tmp_gmt, 1) > 0

            fields = unique([spec(ii).tmp_out_nonsingle.Properties.VariableNames,spec(ii).tmp_out.Properties.VariableNames]);

            for jj = 1:numel(fields)
                try
                    spec(ii).tmp_out_nonsingle.(fields{jj}) = any2str(spec(ii).tmp_out_nonsingle.(fields{jj}));
                catch ME
                    true;
                end

                try
                    spec(ii).tmp_out.(fields{jj}) = any2str(spec(ii).tmp_out.(fields{jj}));
                catch ME
                    true;
                end
            end

            num_moas_represented = num_moas_represented + 1;
        else
            disp(moa_class)
            disp("NOT INCLUDED - NO NON-SINGLETON CLUSTERS")
        end
        
    catch ME
        disp('Error from spectral clustering')
        disp(ME)
    end
    
end

disp('All done')

disp(sprintf('Number of MOAs with non-singleton clusters: %d', num_moas_represented))
disp(sprintf('Total number of MOAs: %d', numel(moas)))

disp(sprintf('Proportion of MOAs with non-singleton clusters: %f', num_moas_represented / numel(moas)))

out_lap = struct2table(out_lap);

out = cat(1, spec.tmp_out);

out_nonsingle = cat(1, spec.tmp_out_nonsingle);

gmt = cat(1, spec.tmp_gmt);

height(out)

headt(out)

height(out_nonsingle)

headt(out_nonsingle)

disp(sprintf('Number of clusters of size greater than 2: %d', sum([gmt.len]>2)))

disp(sprintf('Number of clusters of size 2: %d', sum([gmt.len]<=2)))

disp(sprintf('Number of clusters: %d', size(gmt, 1)))

total_num_trts = numel(unique(out.cid));

num_trts_in_clusters = numel(unique(out_nonsingle.cid));

num_trts_singletons = numel(setdiff(unique(out.cid), unique(out_nonsingle.cid)));

disp(sprintf('Number of treatments: %d', total_num_trts))

disp(sprintf('Number of treatments in clusters: %d', num_trts_in_clusters))

disp(sprintf('Number of singleton, non-cluster treatments: %d', num_trts_singletons))

disp(sprintf('Proportion of treatments in clusters: %f', num_trts_in_clusters / total_num_trts))

disp(sprintf('Proportion of singleton, non-cluster treatments: %f', num_trts_singletons / total_num_trts))

% Save the outputs
wtable(out_lap,fullfile(outdir, moa_k_values_filename))
wtable(out,fullfile(outdir, clusters_tbl_incl_singletons_filename))
wtable(out_nonsingle,fullfile(outdir, clusters_tbl_filename))
mkgmt(fullfile(outdir, clusters_gmt_filename),gmt)

% # LOOCV section

if prepare_loocv

    unique_kabx_cmpds_tbl = rtable(unique_kabx_cmpds_tbl_path);

    size(unique_kabx_cmpds_tbl)
    headt(unique_kabx_cmpds_tbl)
    
    unique_kabx_cmpds_list = unique(unique_kabx_cmpds_tbl.kabx_cmpd);

    length(unique_kabx_cmpds_list)
    
    number_of_cmpds_loocv = length(unique_kabx_cmpds_list)
    
    index_cmpds_loocv = 1:number_of_cmpds_loocv;
    
    if demo_loocv
       if strcmp(demo_loocv_number_or_list, 'number')
           number_of_cmpds_loocv = max(1, demo_loocv_number_cmpds)
           
           index_cmpds_loocv = 1:number_of_cmpds_loocv
           
       elseif strcmp(demo_loocv_number_or_list, 'list')
           number_of_cmpds_loocv = length(demo_loocv_list_cmpds)
           
           index_cmpds_loocv = find(ismember(unique_kabx_cmpds_list, demo_loocv_list_cmpds))'
       else
           error('Invalid input for demo_loocv_number_or_list: number or list')
       end
    end
    
    disp(sprintf('Number of KABX compounds to be processed in LOOCV: %d', length(index_cmpds_loocv)))
    
    for i = index_cmpds_loocv

        % If the current iteration number is a multiple of 50
        if mod(i, 50) == 0
            % Print a status message
            fprintf('Currently at iteration %d\n', i);
        end

        leave_out_cmpd = unique_kabx_cmpds_list(i);
        
        loo_wkdir = fullfile(wkdir, strcat(results_subdir_prefix, strjoin(unique_kabx_cmpds_list(i))));

        mk_cd_dir(loo_wkdir, false);
        
        loo_outdir = fullfile(loo_wkdir, outdir_name)
        
        mk_cd_dir(loo_outdir, false);
        
        % step specific commands
        
        moas = parse_gmt(fullfile(loo_wkdir,moa_gmt_savepath));
        
        spec = {};
        out_lap = {};

        out = {};
        out_nonsingle = {};
        gmt = {};

        num_moas_represented = 0;

        for ii = 1:numel(moas)
            loop_progress(ii,numel(moas),1)

            % Set MoA class
            moa_class = moas(ii).head;
            out_lap(ii).moa_class = moa_class;
            disp(moa_class)

            %[strrep(strrep(moa_class,' ','_'),'/','-'),'_corr_n*.gctx']

            % Load gcts
            g = glob(fullfile(loo_wkdir, moa_gctdir,[strrep(strrep(moa_class,' ','_'),'/','-'),'_corr_n*.gctx']));

            try
                c = parse_gctx(g{1});
            catch ME
                disp(sprintf('Error: %s correlation gctx file not found', moa_class))
                disp(ME)
            end

            g = glob(fullfile(loo_wkdir, moa_gctdir,[strrep(strrep(moa_class,' ','_'),'/','-'),'_corr_rank_n*.gctx']));

            try
                cr = parse_gctx(g{1});
            catch ME
                disp(sprintf('Error: %s correlation rank gctx file not found', moa_class))
                disp(ME)
            end

            moa_size = numel(cr.rid);
            log_moa_size = log(moa_size);

            % Set threshold
            if dynamic_thrsh_per_moa
                thrsh = round(log_moa_size * thrsh_rank);
            else
                thrsh = thrsh_rank*thrsh_factor;
            end

            out_lap(ii).moa_size = moa_size;
            out_lap(ii).log_moa_size = log_moa_size;
            out_lap(ii).thrsh_rank = thrsh_rank;
            out_lap(ii).thrsh_factor = thrsh_factor;
            out_lap(ii).thrsh = thrsh;

            % Run spectral clustering
            try
                [spec(ii).tmp_out_nonsingle,spec(ii).tmp_out,spec(ii).tmp_gmt,L,out_lap(ii).k,en,den,out_lap(ii).k_gap_den,out_lap(ii).k_med_gap_den,out_lap(ii).k_num_zero_plus_one,out_lap(ii).k_num_zero] = spectral_clustering_for_pcls(c,cr,moa_class,thrsh,k_type,loo_outdir,loocv_save_out,loocv_save_fig,rng_seed,show_hclust);
                close all
                
                if size(spec(ii).tmp_gmt, 1) > 0

                    fields = unique([spec(ii).tmp_out_nonsingle.Properties.VariableNames,spec(ii).tmp_out.Properties.VariableNames]);

                    for jj = 1:numel(fields)
                        try
                            spec(ii).tmp_out_nonsingle.(fields{jj}) = any2str(spec(ii).tmp_out_nonsingle.(fields{jj}));
                        catch ME
                            true;
                        end

                        try
                            spec(ii).tmp_out.(fields{jj}) = any2str(spec(ii).tmp_out.(fields{jj}));
                        catch ME
                            true;
                        end
                    end

                    num_moas_represented = num_moas_represented + 1;
                else
                    disp(moa_class)
                    disp("NOT INCLUDED - NO NON-SINGLETON CLUSTERS")
                end
                
            catch ME
                disp('Error from spectral clustering')
                disp(ME)
            end

        end

        disp('All done')

        disp(sprintf('Number of MOAs with non-singleton clusters: %d', num_moas_represented))
        disp(sprintf('Total number of MOAs: %d', numel(moas)))

        disp(sprintf('Proportion of MOAs with non-singleton clusters: %f', num_moas_represented / numel(moas)))
        
        out_lap = struct2table(out_lap);

        out = cat(1, spec.tmp_out);

        out_nonsingle = cat(1, spec.tmp_out_nonsingle);

        gmt = cat(1, spec.tmp_gmt);

        height(out)

        headt(out)

        height(out_nonsingle)

        headt(out_nonsingle)

        disp(sprintf('Number of clusters of size greater than 2: %d', sum([gmt.len]>2)))

        disp(sprintf('Number of clusters of size 2: %d', sum([gmt.len]<=2)))

        disp(sprintf('Number of clusters: %d', size(gmt, 1)))

        total_num_trts = numel(unique(out.cid));

        num_trts_in_clusters = numel(unique(out_nonsingle.cid));

        num_trts_singletons = numel(setdiff(unique(out.cid), unique(out_nonsingle.cid)));

        disp(sprintf('Number of treatments: %d', total_num_trts))

        disp(sprintf('Number of treatments in clusters: %d', num_trts_in_clusters))

        disp(sprintf('Number of singleton, non-cluster treatments: %d', num_trts_singletons))

        disp(sprintf('Proportion of treatments in clusters: %f', num_trts_in_clusters / total_num_trts))

        disp(sprintf('Proportion of singleton, non-cluster treatments: %f', num_trts_singletons / total_num_trts))
        
        
        % Save the outputs
        wtable(out_lap,fullfile(loo_outdir, moa_k_values_filename))
        wtable(out,fullfile(loo_outdir, clusters_tbl_incl_singletons_filename))
        wtable(out_nonsingle,fullfile(loo_outdir, clusters_tbl_filename))
        mkgmt(fullfile(loo_outdir, clusters_gmt_filename),gmt)
        
    
        disp('Done')
        
    end
    
end

disp(sprintf('Number of computational threads: %d', maxNumCompThreads));
maxNumCompThreads('automatic'); % reset number of computational threads to default based on hardware resources
disp(sprintf('Number of computational threads: %d', maxNumCompThreads));