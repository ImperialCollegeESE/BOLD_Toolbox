function [n_optimal_fastchg_given_Tcombo,sim_results_structs_fastchg,param_fastchg_structs, all_combos_Tinit_Treservoir_K, optimal_nomin_cell_capacities_fastchg_vector, exit_reason_flag_vector,ix_saturation_vector,SOC_log_vector] = ...
        compute_optim_layers_fastchg(neg_to_pos_len_ratio,n_optimal_acc,n_max,no_of_cells_in_pack,P_batt_charging_kW, T_init_vector_limits_degC, T_reservoir_limits_degC,...
        temperature_permutations, platform, SOC_init_fast_chg_pct, CP_finish_pct_fast_chg, search_strat)
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

    T_max_charge_Cel = 55; % degC

    if strcmp(platform,'BEV')
        run('common_settings_layer_optimisation_BEV');  % Apply common settings for layer optimisation
    elseif strcmp(platform,'PHEV')
        run('common_settings_layer_optimisation_PHEV'); % Apply common settings for layer optimisation
    end

    T_abs_zero_K         = 273.15;
    sim_datalog_interval = 0.1;

    %% Apply parameters in temperature loop & call function to evaluation number of layers
    P_batt_charging_W = P_batt_charging_kW*10^3;               % Charging power, W
    P_cell            = P_batt_charging_W/no_of_cells_in_pack; % Needs DEBUG!! IDA fails to converge above a certain power limit

    if strcmp(temperature_permutations, 'all')
        T_reservoir_limits_expanded_vector_degC                 = sort([T_init_vector_limits_degC T_reservoir_limits_degC]);    % Combine init & reservoir temps & sort ascending order
        [T_init_vector_degC_meshgrid,T_reservoir_degC_meshgrid] = meshgrid(T_init_vector_limits_degC,T_reservoir_limits_expanded_vector_degC);
        temporary_var = cat(2,T_init_vector_degC_meshgrid',T_reservoir_degC_meshgrid');
        all_combos_Tinit_Treservoir_C = unique(reshape(temporary_var,[],2),'rows'); % Each row in this matrix contains all combinations of initial and reservoir temperatures
    elseif strcmp(temperature_permutations, 'extremes_only')
        middle_temperature = 25; % degC
        all_combos_Tinit_Treservoir_C = zeros(5,2); % Pre-allocate for the four extreme case & the "middle" case. Hence 5 rows
        all_combos_Tinit_Treservoir_C(1,1) = min(T_init_vector_limits_degC);
        all_combos_Tinit_Treservoir_C(2,1) = max(T_init_vector_limits_degC);
        all_combos_Tinit_Treservoir_C(3,1) = middle_temperature;
        all_combos_Tinit_Treservoir_C(4,1) = max(T_init_vector_limits_degC);
        all_combos_Tinit_Treservoir_C(5,1) = min(T_init_vector_limits_degC);

        all_combos_Tinit_Treservoir_C(1,2) = min(T_reservoir_limits_degC);
        all_combos_Tinit_Treservoir_C(2,2) = min(T_reservoir_limits_degC);
        all_combos_Tinit_Treservoir_C(3,2) = middle_temperature;
        all_combos_Tinit_Treservoir_C(4,2) = max(T_reservoir_limits_degC);
        all_combos_Tinit_Treservoir_C(5,2) = max(T_reservoir_limits_degC);
    end

    all_combos_Tinit_Treservoir_K = T_abs_zero_K + all_combos_Tinit_Treservoir_C;
    no_of_temperature_combos      = size(all_combos_Tinit_Treservoir_K,1);

    n_optimal_fastchg_given_Tcombo               = NaN(no_of_temperature_combos,1); % Define an object to hold all of the n_min results for all temperature combinations
    optimal_nomin_cell_capacities_fastchg_vector = NaN(size(n_optimal_fastchg_given_Tcombo));
    sim_results_structs_fastchg                  = cell(no_of_temperature_combos,1);
    param_fastchg_structs                        = cell(no_of_temperature_combos,1);           % Create an object to hold all params used for charging
    exit_reason_flag_vector                      = cell(no_of_temperature_combos,1);
    ix_saturation_vector                         = NaN(size(n_optimal_fastchg_given_Tcombo));
    SOC_log_vector                               = cell(no_of_temperature_combos,1);

    commandwindow;
    if strcmp(search_strat, 'binary')
        cut_point_layers = [];
    end
    %% Determine, for all temperature permutations, the optimal number of layers required for CP fast-charging spec
    for idx_temperature_combo = 1:no_of_temperature_combos
        T_init_K      = all_combos_Tinit_Treservoir_K(idx_temperature_combo,1);
        T_reservoir_K = all_combos_Tinit_Treservoir_K(idx_temperature_combo,2);

        if idx_temperature_combo > 1 && strcmp(search_strat, 'binary')
            if isnan(n_optimal_fastchg_given_Tcombo(idx_temperature_combo-1))
                cut_point_layers = [];
            else
                cut_point_layers = n_optimal_fastchg_given_Tcombo(idx_temperature_combo-1);
            end
        elseif idx_temperature_combo == 1 && strcmp(search_strat, 'linear')
            cut_point_layers = [];                    % Required as a placeholder argument in the call to search_layers_fastchg
        end

        fprintf('\nNow simulating Temperature combo %2d of %2d, with T_init = %4.2f degC & T_reservoir = %4.2f degC\n\n',idx_temperature_combo,no_of_temperature_combos,T_init_K-T_abs_zero_K,T_reservoir_K-T_abs_zero_K);

        param_fastchg{1} = apply_params_for_layer_optim_fcn(SOC_init_fast_chg_pct, ... % init_cell_soc_percent at beginning of fast charge
            OperatingMode, ...                        % OperatingMode - current density input - irrelevant on this call only
            TemperatureEnabled, ...                   % TemperatureEnabled
            cutoffSOC, ...                            % cutoffSOC
            CutoverSOC, ...                           % CutoverSOC
            CutoffVoltage, ...                        % CutoffVoltage
            CutoverVoltage, ...                       % CutoverVoltage
            solidDiffusionMode, ...                   % solidDiffusionMode
            (T_init_K), ...                           % T_init - updated for all permutations
            (T_reservoir_K), ...                      % T_reservoir - updated for all permutations
            (T_abs_zero_K + T_max_charge_Cel), ...    % Tmax (kelvin)
            sim_datalog_interval, ...                 % sim_datalog_interval
            axial_nodes, ...                          % axial_nodes  (i.e. in the through-thickness direction)
            radial_nodes);                            % radial_nodes (i.e. no. of shells for discretising a spherical particle)

        param_fastchg{1}.enable_csneg_Saturation_limit = 1; % This parameter, when set to 1, enforces simulation termination when surface concentration of any node in the neg electrode hits the saturation value (set in parameter above)
        param_fastchg{1}.cs_neg_saturation = ((0.01*param_fastchg{1}.CutoverSOC*(param_fastchg{1}.theta_max_neg-param_fastchg{1}.theta_min_neg) + param_fastchg{1}.theta_min_neg))*param_fastchg{1}.cs_maxn;

        [n_optimal_fastchg_given_Tcombo(idx_temperature_combo), sim_results_structs_fastchg{idx_temperature_combo},param_fastchg_structs{idx_temperature_combo},optimal_nomin_cell_capacities_fastchg_vector(idx_temperature_combo), exit_reason_flag_vector{idx_temperature_combo}, ix_saturation_vector(idx_temperature_combo), SOC_log_vector{idx_temperature_combo}] ...
            = search_layers_fastchg(n_optimal_acc,n_max,cut_point_layers,neg_to_pos_len_ratio,P_cell,CP_finish_pct_fast_chg,param_fastchg, search_strat);
        fprintf('Completed optimal layer choice computation for fast-charging with T_init = %4.2f degC & T_reservoir = %4.2f degC\n\n',T_init_K-T_abs_zero_K, T_reservoir_K-T_abs_zero_K);
        fprintf('\nFor all temperature combos up to & including this one (%2d of %2d), the optimal layer choice for fast-charging is % d\n',idx_temperature_combo,no_of_temperature_combos,max(n_optimal_fastchg_given_Tcombo(1:idx_temperature_combo)));
        fprintf('---------------------------------------------------------------------------------------------------------------------\n\n');
    end     % End of temperature-stepping loop
    % Fast-charging upto 80% SOC (spec) is performed and optimal number of layers determined. Below code is optional.

    %% The user may choose to uncomment the following lines of code, if fast-charging upto 100% is desired. This does not influence the optimal layer choice computed above
    % bang_bang_threshold_pct = 99;   % in percent, currently fixed to this value for all charging powers.

    % Minimum num_layers has been determined for constant power to saturation. Complete pulsed power charging phase with this num_layers
    %
    % % Inform user of status of stage 4, pulsed testing
    % P_density = P_cell/param_fastchg{1}.overall_surface_area_for_given_layers;  % Re-compute P_density using surface area for highest number of layers
    % param_100pct_soc = Parameters_init_suppliedSOC_pct(100); % optional requirement (only for compute_capacity_for_layer_fcn)
    % param_0pct_soc = Parameters_init_suppliedSOC_pct(0);     % optional requirement (only for compute_capacity_for_layer_fcn)
    % neg_electrode_capacity_Ah = compute_capacity_for_layer_fcn(n_optimal_fastchg, param_100pct_soc, param_0pct_soc, param_fastchg{1}.len_n);
    % cell_capacity_Ah = neg_electrode_capacity_Ah;
    % fprintf('|   Pulse power  |')
    % fprintf('   %d   |       %.2f       |       %.2f      |       %.2f        | In progress...\n\n', n_optimal_fastchg, param_fastchg{1}.overall_surface_area_for_given_layers, P_density, cell_capacity_Ah)
    % param_fastchg{1}.suppress_status_prints = 0;    % Reset LIONSIMBA status updates to allow output to command window
    % param_fastchg{1}.sim_datalog_interval = 1;   % Required
    %
    % % Obtain end state of stage 3 to continue charging from
    % % temp_results_CP_charge_end = results_fastchg(end);
    % t_local_start_pulsing = results_fastchg(end).time{1}(end);   % To continue from the constant power charge of stage 3 above
    % CP_charge_time_pulsing = results_fastchg(end).time{1}(end);  % Save for displaying at end of simulation
    % sampling_interval_pulsing = 1;                                    % timestep
    % initialState_pulsing = results_fastchg(end).initialState;    % Update the initial states
    % t_local_finish_pulsing = t_local_start_pulsing + sampling_interval_pulsing;
    %
    % max_cs_surface_neg_electrode_pulsing = max(results_fastchg(end).cs_surface{1}(end,param_fastchg{1}.Np+1:end));
    % cell_soc_pct = results_fastchg(end).SOC{1}(end);
    % bang_bang_lower_cs_limit = (bang_bang_threshold_pct/100)*cs_neg_saturation;
    % hit_lower_threshold_flag = 0;  % flag needed to implement hysterisis (bang-bang control)
    % param_fastchg{1}.RelTol = 1e-6;
    % param_fastchg{1}.AbsTol = 1e-6;
    %
    % % Pulse intermittently until upper SOC cutover is reached
    % while cell_soc_pct <= param_fastchg{1}.CutoverSOC
    %     if hit_lower_threshold_flag==0 && max_cs_surface_neg_electrode_pulsing>=cs_neg_saturation
    %         P_density_local = 1e-1;     % zero power density causes IDA to fail. Approximate zero with 1e-2 Watts
    %         hit_lower_threshold_flag = 1;
    %     elseif hit_lower_threshold_flag==1 && max_cs_surface_neg_electrode_pulsing<=bang_bang_lower_cs_limit
    %         P_density_local = P_density;
    %         hit_lower_threshold_flag = 0;
    %     elseif hit_lower_threshold_flag==0 && max_cs_surface_neg_electrode_pulsing<=cs_neg_saturation
    %         P_density_local = P_density; % P_density_local = P_density_local;
    %     else
    %
    %     end
    %
    %     try
    %         results_local_pulsing = startSimulation(t_local_start_pulsing,t_local_finish_pulsing,initialState_pulsing,P_density_local,param_fastchg);
    %     catch ME
    %         fprintf('Encountered an error during pulse charging\n');
    %         fprintf(ME.identifier); fprintf('\n');
    %     end
    %
    %     max_cs_surface_neg_electrode_pulsing = max(results_local_pulsing.cs_surface{1}(end,param_fastchg{1}.Np+1:end));
    %     t_local_start_pulsing = results_local_pulsing.time{1}(end);       % Update the starting integration time;                         % start_time proceeds forward by the sampling interval of the drive cycle
    %     t_local_finish_pulsing = t_local_start_pulsing + sampling_interval_pulsing;     % start_time proceeds forward by the sampling interval of the drive cycle
    %     cell_soc_pct = results_local_pulsing.SOC{1}(end);
    %     initialState_pulsing = results_local_pulsing.initialState;  % Update the initial states
    %     fprintf('Total charge time: %4.1f s | SOC: %4.2f %% | Target SOC: %4.2f %% | Surf. lithiation: %4.2f %% saturated | Power density: %3.2f | Temp: %5.2f degC |\n', t_local_finish_pulsing, results_local_pulsing.SOC{1}(end), param_fastchg{1}.CutoverSOC, (max_cs_surface_neg_electrode_pulsing/cs_neg_saturation)*100, P_density_local, max(results_local_pulsing.Temperature{1}(end))-273.15);
    %
    %     if results_local_pulsing.exit_reason ~= 0 && results_local_pulsing.exit_reason ~= 4 % break out of while loop & start next iteration of for loop, because the simulation didn't pass this stage
    %         phase_four_result = 'Fail';     % At this point, SOC target wasn't reached because simulation ended unsuccessfully
    %         error = exit_code_to_string(results_local_pulsing.exit_reason);
    %         fprintf('|   Pulse power  |   %d   |       %.2f       |       %.2f      |       %.2f        | %s | %s |\n\n', num_layers, param_fastchg{1}.overall_surface_area_for_given_layers, P_density, cell_capacity_Ah, phase_four_result, error);                      %
    %         break;  % Exit the while loop, because this num_layers has failed for reason: results_local.exit_reason
    %     end
    % end     % End of SOC-while loop
    %
    % if cell_soc_pct >= param_fastchg{1}.CutoverSOC  % If above loop exited "naturally", this condition should be satisfied
    %     phase_four_result = 'Pass';                             % At this point, "while" loop ended because target SOC was reached
    %     error = exit_code_to_string(results_local_pulsing.exit_reason);
    %     fprintf('|   Pulse power  |   %d   |       %.2f       |       %.2f      |       %.2f        | %s | %s |\n\n', n_optimal_fastchg, param_fastchg{1}.overall_surface_area_for_given_layers, P_density, cell_capacity_Ah, phase_four_result, error);
    %     fprintf('Constant power charge time was: %.2f', CP_charge_time_pulsing);
    %     fprintf('Pulse power charge time was: %.2f', (t_local_finish_pulsing - CP_charge_time_pulsing));
    %     fprintf('Total charge time was: %.2f', t_local_finish_pulsing);
    %     diary('Charge_log');
    % else
    %     fprintf('There was an error with pulse charging')
    % end
    %
    % fprintf('\nFinished charging.');
end     % End of function

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: