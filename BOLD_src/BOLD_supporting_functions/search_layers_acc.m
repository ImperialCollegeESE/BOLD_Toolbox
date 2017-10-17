%% Perform binary search to determine the optimal number of layers for acceleration spec
function [optimal_layers_acc, optimal_results_acc_struct,optimal_param_acc_struct,nominal_cell_capacity,P_wheels_cruise, exit_reason_flag_vector] = ...
        search_layers_acc(n_min, n_max, cut_point_layers, neg_to_pos_len_ratio,Mv_with_pack_overhead,Cd,Av,tf,v_base,eta_drivetrain,no_of_cells_in_pack,vf,param_acc, search_strat)
    % Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
    %                : Davide M. Raimondo, University of Pavia
    % copyright year : 2017
    % Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
    % Licensed       : MIT License

    if strcmp(search_strat, 'binary')
        max_number_of_iterations      = ceil(log2(n_max - n_min  + 2));      % worst case iterations for binary search
        sim_results_structs_acc       = cell(max_number_of_iterations,1);
        param_structs_acc             = cell(max_number_of_iterations,1);
        layer_eval_result_flag_vector = zeros(max_number_of_iterations,1);
        nominal_cell_capacity_vector  = zeros(max_number_of_iterations,1);
        P_wheels_cruise_vector        = zeros(max_number_of_iterations,1);
        exit_reason_flag_vector       = cell(max_number_of_iterations,1);

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

            fprintf('\nBinary search iteration no. %d begins: Attempting to begin acceleration run with %d layers ...\n',search_iter,num_layer_trial);
            [layer_eval_result_flag_vector(search_iter),sim_results_structs_acc{search_iter},param_structs_acc{search_iter},nominal_cell_capacity_vector(search_iter),P_wheels_cruise_vector(search_iter)] = ...
                eval_layer_feasibility(num_layer_trial, neg_to_pos_len_ratio,Mv_with_pack_overhead,Cd,Av,tf,v_base,eta_drivetrain,no_of_cells_in_pack,vf,param_acc);

            if isempty(sim_results_structs_acc{search_iter})
                exit_reason_flag_vector{search_iter} = NaN;
            else
                exit_reason_flag_vector{search_iter} = sim_results_structs_acc{search_iter}.exit_reason;
            end

            layer_eval_result_flag = layer_eval_result_flag_vector(search_iter);
            if (layer_eval_result_flag == 0) && num_layer_trial ~= (n_max-1)          % num_layer_trial failed, but we haven't yet reached the upper bound (R) on the search space
                L = m;                                                                % Increment the lower bound
                optimal_layers_acc = num_layer_trial + 1;                             % This can only be incremented if the choice m isn't already == n_max
            elseif (layer_eval_result_flag == 1)
                R = m;
                optimal_layers_acc = num_layer_trial;
            elseif (layer_eval_result_flag == 0) && num_layer_trial == (n_max-1)      % num_layer_trial failed & we HAVE reached the upper bound (R) on the search space. i.e. if this condition passes, there was no suitable number of layers in the range n_min to n_max for the temperature combination
                L = m;                                                                % Increment the lower bound - necessary here or while loop will be infinite
                optimal_layers_acc = NaN;                                             % No suitable layer number in the range was successful
            else
                fprintf('Something unexpected happened!');
            end
        end
        % layer_eval_result_flag_vector(isnan(layer_eval_result_flag_vector))=[];
        optimal_iter = find(layer_eval_result_flag_vector,1,'last');                % Find the last 1 element that is non-zero in layer_eval_result_flag_vector

        if isempty(optimal_iter)                    % find will return an empty array if none of the layer choice were successful
            fprintf('Warning! layer_eval_result_flag_vector contains only zeros - none of the attempted layers were successful!\n');
            optimal_results_acc_struct = NaN;
            optimal_param_acc_struct   = NaN;
            nominal_cell_capacity      = NaN;
            P_wheels_cruise            = NaN;
        else                                        % If one of the layers choices was a success, return the results & param from that call to startSimulation
            optimal_results_acc_struct = sim_results_structs_acc(optimal_iter);         % Use the above as the index
            optimal_param_acc_struct   = param_structs_acc(optimal_iter);
            nominal_cell_capacity      = nominal_cell_capacity_vector(optimal_iter);
            P_wheels_cruise            = P_wheels_cruise_vector(optimal_iter);
        end

    else    % Linear search through layers
        search_iter                   = 1;
        num_layer_trial               = n_min;
        sim_results_structs_acc       = cell(length(n_min:n_max),1);          % Over-allocation, since don't know when "while" loop will end
        param_structs_acc             = cell(length(n_min:n_max),1);                % Over-allocation, since don't know when "while" loop will end
        nominal_cell_capacity_vector  = zeros(length(n_min:n_max),1);
        layer_eval_result_flag_vector = zeros(length(n_min:n_max),1);   % Over-allocation, since don't know when "while" loop will end
        P_wheels_cruise_vector        = zeros(length(n_min:n_max),1);

        exit_reason_flag_vector = cell(length(n_min:n_max),1);         % Over-allocation, since don't know when "while" loop will end

        while num_layer_trial <= n_max && layer_eval_result_flag_vector(search_iter) ~= 1   % Ends when n_max is reached, or when layer is a success
            fprintf('\nLinear search iteration no. %d begins: Attempting to begin acceleration run with %d layers...\n',search_iter,num_layer_trial);

            if num_layer_trial > n_min          % Prevents search_iter getting incremented on the first pass
                search_iter = search_iter + 1;
            end

            [layer_eval_result_flag_vector(search_iter),sim_results_structs_acc{search_iter},param_structs_acc{search_iter},nominal_cell_capacity_vector(search_iter),P_wheels_cruise_vector(search_iter)] = ...
                eval_layer_feasibility(num_layer_trial, neg_to_pos_len_ratio,Mv_with_pack_overhead,Cd,Av,tf,v_base,eta_drivetrain,no_of_cells_in_pack,vf,param_acc);

            if isempty(sim_results_structs_acc{search_iter})
                exit_reason_flag_vector{search_iter} = NaN;
            else
                exit_reason_flag_vector{search_iter} = sim_results_structs_acc{search_iter}.exit_reason;
            end

            fprintf('layer_eval_result_flag_vector value: %d, exit_reason: %d \n\n\n', layer_eval_result_flag_vector(search_iter), exit_reason_flag_vector{search_iter});

            num_layer_trial = num_layer_trial + 1;

        end % End of "for" loop in linear layer search

        if layer_eval_result_flag_vector(search_iter) == 1  % Above while loop can exit because n_max is reached. Need to check if that occurred
            optimal_layers_acc         = num_layer_trial - 1;
            optimal_results_acc_struct = sim_results_structs_acc(search_iter);
            optimal_param_acc_struct   = param_structs_acc(search_iter);
            nominal_cell_capacity      = nominal_cell_capacity_vector(search_iter);
            P_wheels_cruise            = P_wheels_cruise_vector(search_iter);
        else
            fprintf('Warning! layer_eval_result_flag_vector contains only zeros - none of the attempted layers were successful!\n')
            optimal_layers_acc         = NaN;
            optimal_results_acc_struct = NaN;
            optimal_param_acc_struct   = NaN;
            nominal_cell_capacity      = NaN;
            P_wheels_cruise            = NaN;
        end

    end     % End of "if" for search_strat search selected

end % End of 'binary_search_layers_acc_simplified' function

function [layer_eval_result_flag,results_acc_constPower,param_acc_simulation,nominal_cell_capacity,P_wheels_cruise] = eval_layer_feasibility(num_layers, neg_to_pos_len_ratio,Mv_with_pack_overhead,Cd,Av,tf_acc,v_base,eta_drivetrain,no_of_cells_in_pack,vf,param_acc_original)
    param_acc_simulation = param_acc_original;
    results_acc_constPower  = [];
    % nominal_cell_capacity = 'not computed' ; % for optionally printing cell capacity later

    %% Re-compute parameters which are functions of number of layers
    [~, len_p, len_n] = compute_domain_thicknesses(num_layers, neg_to_pos_len_ratio, param_acc_original);
    param_acc_simulation{1}.len_p = len_p;
    param_acc_simulation{1}.len_n = len_n;

    param_acc_simulation{1}.overall_surface_area_for_given_layers = num_layers*param_acc_original{1}.surface_area_per_face_Northrop_cell;
    [param_acc_simulation{1},~,~,~] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(num_layers, param_acc_simulation{1}); % cell mass and Cp are appended to param struct

    Mv                 = Mv_with_pack_overhead + param_acc_simulation{1}.mass_cell*no_of_cells_in_pack;  % Entire vehicle mass with pack (neglecting rot. components)
    P_wheels_cruise    = compute_cruise_power(vf, Mv, Cd, Av);    % Vehicle's power demand (watts) without considering power required to accelerate its mass.
    P_acc_mass         = 0.5*Mv*(v_base^2 + vf^2)/tf_acc;        % Power required solely for accelerating the vehicle's mass
    P_wheels_acc_total = P_wheels_cruise + P_acc_mass;    % Total power demand during acceleration phase (0-60mph)
    % disp(P_wheels_acc_total/1000);
    % fprintf('\n');

    P_batt_acc    = -(P_wheels_acc_total/eta_drivetrain);  % Battery pack power (negative sign introduced here for adhering to LIONSIMBA conventions, i.e. discharge power is negative. Accelerating discharges the pack)
    P_cell_acc    = P_batt_acc/no_of_cells_in_pack;        % Power to be supplied by a single cell
    P_density_acc = P_cell_acc/param_acc_simulation{1}.overall_surface_area_for_given_layers;

    fprintf('\nWith %3d layer(s), acceleration power is %4.2f W per cell. This results in applied power density of %6.2f W/m^2 for simulating this layer choice.\n',num_layers,P_cell_acc,P_density_acc);
    fprintf('P_batt_acc = %f kW\n', P_batt_acc/1000);
    %% Optionally compute cell capacity for this layer choice
    param_100pct_soc      = Parameters_init_suppliedSOC_pct(100);
    param_0pct_soc        = Parameters_init_suppliedSOC_pct(0);
    nominal_cell_capacity = compute_capacity_for_layer_fcn(num_layers,param_100pct_soc,param_0pct_soc,len_n);

    %% Attempt constant power acceleration (i.e. discharge) simulation for 'tf' seconds
    try
        fprintf('\nConst. power acceleration simulation in progress ........\n');
        results_acc_constPower = startSimulation(0,tf_acc,[],P_density_acc,param_acc_simulation);
        fprintf('\nTermination Voltage: %4.3f | Cell Capacity: %4.2f Ah | Cell SOC: %5.2f %% | T_init: %4.2f degC | T_ambient: %4.2f degC | Cell Temp.: %4.2f degC\n', results_acc_constPower.Voltage{1}(end),nominal_cell_capacity,results_acc_constPower.SOC{1}(end),param_acc_original{1}.T_init-273.15,param_acc_original{1}.Tref-273.15,max(results_acc_constPower.Temperature{1}(end))-273.15);

        if ismember(1,results_acc_constPower.exit_reason) || ismember(3,results_acc_constPower.exit_reason) || ismember(5,results_acc_constPower.exit_reason) % hit cutoff voltage, cutoff SOC or upper temperature limit
            layer_eval_result_flag = 0;

            if ismember(1,results_acc_constPower.exit_reason)
                fprintf('\nHit Cutoff Voltage! ...\n');
            end

            if ismember(3,results_acc_constPower.exit_reason)
                fprintf('\nHit Cutoff SOC! ...\n');
            end

            if ismember(5,results_acc_constPower.exit_reason)
                fprintf('\nExceeded safe operating cell temperature during acceleration run ...\n');
            end

            fprintf('\nHence unsuccessful in meeting acceleration spec with %d layers ...\n\n\n',num_layers);
        else
            fprintf('\nSuccess with %3d layers! Completed acceleration run without hitting cutoff voltage, cutoff SOC or temperature limit...\n\n',num_layers);
            layer_eval_result_flag = 1;
        end

    catch
        fprintf('\nToo high power density for IDA (i.e. the solver) to converge\n');
        fprintf('\nUnsuccessful with %2d number of layers ...\n\n\n',num_layers);
        layer_eval_result_flag = 0;
    end
end

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: