%% Perform search to determine the optimal number of layers for fast-charging
function [optimal_layers_fastchg, optimal_results_fastchg_struct,optimal_param_fastchg_struct,nominal_cell_capacity,exit_reason_flag_vector, ix_saturation, SOC_log_vector] = ...
    search_layers_fastchg(n_min, n_max, cut_point_layers, neg_to_pos_len_ratio,P_cell,CP_finish_pct_fast_chg,param_fastchg, search_strat)
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

% Regarding binary search option only: The key assumption of monotonicity of result vector needs
% to be maintained. There is an important consideration here.  layer_eval_result_flag_vector is
% assumed as a proxy for whether a given layer is "successful or not". This assumes that the
% applied power density is high enough to saturate the outer layer, but not exceed the desired
% fast-charge SOC (i.e. the codes have a subtle, inherent bias towards low layer counts)
% Note that even if cutOver SOC is exceeded, we deem the boolean 'layer_eval_result_flag'
% to be 1, solely due to the above reasoning.

if strcmp(search_strat, 'binary')
    max_number_of_iterations      = ceil(log2(n_max - n_min  + 2));    % worst case iterations for binary search
    sim_results_structs_fastchg   = cell(max_number_of_iterations,1);
    param_structs_fastchg         = cell(max_number_of_iterations,1);
    layer_eval_result_flag_vector = zeros(max_number_of_iterations,1);
    nominal_cell_capacity_vector  = zeros(max_number_of_iterations,1);
    exit_reason_flag_vector       = cell(max_number_of_iterations,1);
    ix_saturation_vector          = zeros(max_number_of_iterations,1);
    SOC_log_vector                = zeros(max_number_of_iterations,1);

    search_iter = 0;
    L = 1;
    R = n_max - n_min + 1;
    while (L < R-1)
        search_iter = search_iter + 1;
        if search_iter == 1 && ~isempty(cut_point_layers)
            m = cut_point_layers - n_min + 1;
        else
            m = round((L + R)/2);
        end
        num_layer_trial = m + n_min - 1;

        fprintf('\nBinary search iteration no. %d begins: Attempting to begin fast-charging with %d layers ...\n',search_iter,num_layer_trial);
        [layer_eval_result_flag_vector(search_iter),sim_results_structs_fastchg{search_iter},param_structs_fastchg{search_iter},nominal_cell_capacity_vector(search_iter),ix_saturation_vector(search_iter)] = ...
            eval_layer_feasibility(num_layer_trial, neg_to_pos_len_ratio, P_cell, CP_finish_pct_fast_chg, param_fastchg);

        if isempty(sim_results_structs_fastchg{search_iter})    % Log the exit reason, if it exists
            exit_reason_flag_vector{search_iter} = NaN;
        else
            exit_reason_flag_vector{search_iter} = sim_results_structs_fastchg{search_iter}.exit_reason;
            SOC_log_vector(search_iter) = sim_results_structs_fastchg{search_iter}.SOC{1}(ix_saturation_vector(search_iter));   % Record the SOC for that layer count when the anode surface was saturated
        end

        layer_eval_result_flag = layer_eval_result_flag_vector(search_iter);
        if (layer_eval_result_flag == 0)
            L = m;
            optimal_layers_fastchg = num_layer_trial + 1;
        elseif (layer_eval_result_flag == 1)
            R = m;
            optimal_layers_fastchg = num_layer_trial;
        end
    end
    optimal_iter = find(layer_eval_result_flag_vector,1,'last');

    if isempty(optimal_iter)                    % find will return an empty array if none of the layer choice were successful
        fprintf('Warning! layer_eval_result_flag_vector contains only zeros - none of the attempted layers were successful!\n');
        optimal_results_fastchg_struct = NaN;
        optimal_param_fastchg_struct = NaN;
        nominal_cell_capacity = NaN;
        ix_saturation = NaN;
    else                                        % If one of the layers choices was a success, return the results & param from that call to startSimulation
        optimal_results_fastchg_struct = sim_results_structs_fastchg(optimal_iter);
        optimal_param_fastchg_struct = param_structs_fastchg(optimal_iter);
        nominal_cell_capacity = nominal_cell_capacity_vector(optimal_iter);
        ix_saturation = ix_saturation_vector(search_iter);                  % Use instead of "end" when indexing the results from time t0 to time "end"
    end

else    % Linear search through layers
    search_iter                   = 1;
    num_layer_trial               = n_min;
    sim_results_structs_fastchg   = cell(length(n_min:n_max),1);
    param_structs_fastchg         = cell(length(n_min:n_max),1);
    layer_eval_result_flag_vector = zeros(length(n_min:n_max),1);
    nominal_cell_capacity_vector  = zeros(length(n_min:n_max),1);
    exit_reason_flag_vector       = cell(length(n_min:n_max),1);
    ix_saturation_vector          = zeros(length(n_min:n_max),1);
    SOC_log_vector                = zeros(length(n_min:n_max),1);

    while num_layer_trial <= n_max && layer_eval_result_flag_vector(search_iter) ~= 1   % Ends when n_max is reached, or when layer is a success

        if num_layer_trial > n_min          % Prevents search_iter getting incremented on the first pass
            search_iter = search_iter + 1;
        end
        fprintf('\nLinear search iteration no. %d begins: Attempting to begin fast-charging with %d layers...\n',search_iter,num_layer_trial);

        [layer_eval_result_flag_vector(search_iter),sim_results_structs_fastchg{search_iter},param_structs_fastchg{search_iter},nominal_cell_capacity_vector(search_iter), ix_saturation_vector(search_iter)] = ...
            eval_layer_feasibility(num_layer_trial, neg_to_pos_len_ratio, P_cell, CP_finish_pct_fast_chg, param_fastchg);

        if isempty(sim_results_structs_fastchg{search_iter})        % Log the exit reason, if it exists
            exit_reason_flag_vector{search_iter} = NaN;
        else
            exit_reason_flag_vector{search_iter} = sim_results_structs_fastchg{search_iter}.exit_reason;
            SOC_log_vector(search_iter) = sim_results_structs_fastchg{search_iter}.SOC{1}(ix_saturation_vector(search_iter));   % Record the SOC for that layer count when the anode surface was saturated
        end

        fprintf('layer_eval_result_flag_vector value: %d, exit_reason: %d \n\n\n', layer_eval_result_flag_vector(search_iter), exit_reason_flag_vector{search_iter});

        num_layer_trial = num_layer_trial + 1;

    end % End of "for" loop in linear layer search

    if layer_eval_result_flag_vector(search_iter) == 1  % Above while loop can exit because n_max is reached. Need to check if that occurred
        optimal_layers_fastchg         = num_layer_trial - 1;
        optimal_results_fastchg_struct = sim_results_structs_fastchg(search_iter);
        optimal_param_fastchg_struct   = param_structs_fastchg(search_iter);
        nominal_cell_capacity          = nominal_cell_capacity_vector(search_iter);
        ix_saturation                  = ix_saturation_vector(search_iter);                  % Use instead of "end" when indexing the results from time t0 to time "end"
    else
        fprintf('Warning! layer_eval_result_flag_vector contains only zeros - none of the attempted layers were successful!\n')
        optimal_layers_fastchg         = NaN;
        optimal_results_fastchg_struct = NaN;
        optimal_param_fastchg_struct   = NaN;
        nominal_cell_capacity          = NaN;
        ix_saturation                  = NaN;
    end

end     % End of "if" for search_strat search selected

end % End of 'binary_search_layers_fastchg_simplified' function

function [layer_eval_result_flag,results_fastchg_constPower,param_fastchg_simulation,nominal_cell_capacity,ix_saturation] = eval_layer_feasibility(num_layers, neg_to_pos_len_ratio, P_cell, CP_finish_pct_fast_chg, param_fastchg_original)
param_fastchg_simulation   = param_fastchg_original;
runtime_local              = 3600;   % seconds
results_fastchg_constPower = [];
% nominal_cell_capacity = 'not computed' ; % for optionally printing cell capacity later

%% Re-compute parameters which are functions of number of layers
[~, len_p, len_n] = compute_domain_thicknesses(num_layers, neg_to_pos_len_ratio, param_fastchg_original);
param_fastchg_simulation{1}.len_p = len_p;
param_fastchg_simulation{1}.len_n = len_n;
param_fastchg_simulation{1}.overall_surface_area_for_given_layers = num_layers*param_fastchg_original{1}.surface_area_per_face_Northrop_cell;
[param_fastchg_simulation{1},~,~,~] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(num_layers, param_fastchg_simulation{1}); % cell mass and Cp are appended to param struct

P_density = P_cell/param_fastchg_simulation{1}.overall_surface_area_for_given_layers;
fprintf('\nWith %3d layer(s), charging power of %4.2f W per cell results in applied power density of %6.2f W/m^2 for simulating this layer choice.\n',num_layers,P_cell,P_density);

%% Compute nominal cell capacity for this layer choice
param_100pct_soc = Parameters_init_suppliedSOC_pct(100);
param_0pct_soc = Parameters_init_suppliedSOC_pct(0);
nominal_cell_capacity = compute_capacity_for_layer_fcn(num_layers,param_100pct_soc,param_0pct_soc,len_n);

%% Attempt constant power simulation until surface concentration exceeds the safe saturation limit

try
%     param_fastchg_simulation{1}.Scope = 1;    % Option for debugging only
    fprintf('\nConst. power charging simulation in progress ........\n');
    results_fastchg_constPower = startSimulation(0,runtime_local,[],P_density,param_fastchg_simulation);

    ix_saturation = find(results_fastchg_constPower.cs_surface{1}(:,param_fastchg_simulation{1}.Np+1:param_fastchg_simulation{1}.Np+1)...
        > param_fastchg_original{1}.cs_sat_thresh*param_fastchg_original{1}.cs_neg_saturation,1); % Determine the vector index where outer shell became saturated - accounts for solver inertia
    if isempty(ix_saturation)   % If the sim didn't run so far that cs_surf reached saturation, no need to trim. Hence, just take the values at the end of the simulation
        ix_saturation = length(results_fastchg_constPower.time{1});    % Equivalent of using "end" in place of ix_saturation when indexing the other variables from results
    end

    max_cs_surface_neg = max(results_fastchg_constPower.cs_surface{1}(ix_saturation,param_fastchg_simulation{1}.Np+1:end)); % Obtain concentration at the time (row) where outer shell became saturated
    fprintf('Cell Capacity: %4.2f Ah | Cell SOC: %5.2f %% | Surf. Conc.: %6.2f %% of safe saturation limit | T_init: %4.2f degC | T_reservoir: %4.2f degC | Cell Temp.: %4.2f degC\n', nominal_cell_capacity,results_fastchg_constPower.SOC{1}(ix_saturation),max_cs_surface_neg*100/param_fastchg_original{1}.cs_neg_saturation,param_fastchg_original{1}.T_init-273.15,param_fastchg_original{1}.Tref-273.15,max(results_fastchg_constPower.Temperature{1}(ix_saturation))-273.15);

    if results_fastchg_constPower.SOC{1}(ix_saturation) < CP_finish_pct_fast_chg || ismember(2,results_fastchg_constPower.exit_reason) || ismember(5,results_fastchg_constPower.exit_reason)
        layer_eval_result_flag = 0;

        if results_fastchg_constPower.SOC{1}(ix_saturation) < CP_finish_pct_fast_chg
            fprintf('\nUnable to attain %5.2f %% SOC (target)\n',CP_finish_pct_fast_chg);
        end

        if ismember(2,results_fastchg_constPower.exit_reason)
            fprintf('\nExceeded cutOver voltage in trying to fast-charge ...\n');
        end

        if ismember(5,results_fastchg_constPower.exit_reason)
            fprintf('\nExceeded safe operating cell temperature in trying to fast-charge ...\n');
        end

        fprintf('\nHence unsuccessful in fast-charging with %d layers ...\n\n\n',num_layers);
    else
        fprintf('\nSuccess with %3d layers! \nReached (or exceeded) fast-charging SOC target of %4.2f %% with %d layers without hitting upper cutoff voltage or temperature ...\n\n\n',num_layers,CP_finish_pct_fast_chg,num_layers);
        layer_eval_result_flag = 1;
    end

catch
    fprintf('\nToo high power density for IDA (i.e. the solver) to converge\n');
    fprintf('\nUnsuccessful with %2d number of layers\n',num_layers);
    layer_eval_result_flag = 0;
    ix_saturation = NaN;        % Required since it's being returned from the function
end
end

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: