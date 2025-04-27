function main(len_qu, rand_seed, c_num)
    % Convert input arguments to numeric values
    l = str2double(len_qu);
    rs = str2double(rand_seed);
    copy_num = str2double(c_num);
    
    % Define data paths
    data_path = fullfile('data', ['l', num2str(l)]);
    rs_str = num2str(rs);
    data_name = fullfile(data_path, [rs_str, '.csv']);
    output_dir = fullfile(data_path, rs_str);
    
    % Create the output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Initialize parameters
    p_num = 64;
    s_ep_list_all = zeros(copy_num, p_num, 'double'); % Preallocate for efficiency
    
    % Determine the optimal number of workers (avoid over-allocation)
    num_cores = feature('numcores');
    num_workers = min(copy_num, num_cores);
    
    % Configure and initialize the parallel pool
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= num_workers
        if ~isempty(pool)
            delete(pool); % Close existing pool if it doesn't match the desired size
        end
        cluster = parcluster('local');
        cluster.JobStorageLocation = strcat('./',data_path,'/',num2str(rs));
    	cluster.NumThreads = 5; % Matches --cpus-per-task=6
        cluster.NumWorkers = num_workers;
        parpool(cluster);
    end

    p_values = (1:p_num) /p_num;
    
    % Parallel loop over the number of copies
    parfor copy_idx = 1:copy_num
        rng(copy_idx + rs * copy_num, 'twister');
        s_ep_list_all(copy_idx, :) = arrayfun(@(pn) rk_diff(l, pn), p_values);
    end
    mean_s_ep = mean(s_ep_list_all, 1);
    writematrix(mean_s_ep, data_name);
    
    % Clean up the parallel pool if it was created within this function
    pool = gcp('nocreate');
    if ~isempty(pool) && pool.NumWorkers == num_workers
        delete(pool);
    end
end
