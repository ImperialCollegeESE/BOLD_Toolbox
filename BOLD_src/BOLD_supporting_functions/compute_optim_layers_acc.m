function [n_optimal_acc_given_Tcombo,sim_results_structs_acc,param_acc_structs, all_combos_Tinit_Treservoir_K_acc, optimal_nomin_cell_capacities_acc_vector,P_wheels_cruise_vector, exit_reason_flag_vector] = ...
        compute_optim_layers_acc(neg_to_pos_len_ratio,n_min,n_max,Mv_with_pack_overhead,Cd,Av,v_base,eta_drivetrain,no_of_cells_in_pack,...
        T_init_vector_limits_degC, T_reservoir_limits_degC, platform, vf, tf, search_strat)
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

    SOC_init_acc_pct = 40;    % Worst case starting condition for acceleration tests (as per SAE spec J1666)

    if strcmp(platform,'BEV')
        run('common_settings_layer_optimisation_BEV');  % Apply common settings for layer optimisation
    elseif strcmp(platform,'PHEV')
        run('common_settings_layer_optimisation_PHEV'); % Apply common settings for layer optimisation
    end

    sim_datalog_interval = 10e-3;   % Choose a higher-data data-logging rate than default, if results.time is used as a stopping criterion

    Tmax_degC                = 55;  % ref: "A review on the key issues for lithium-ion battery management in electric vehicles", Languang Lu, Xuebing Han, Jianqiu Li, Jianfeng Hua, Minggao Ouyang, Journal of Power Sources, 226, 2013, 272 - 288
    T_abs_zero_K             = 273.15;
    T_max_K                  = Tmax_degC + T_abs_zero_K;
    temperature_permutations = 'extremes_only';

    if strcmp(temperature_permutations, 'all')
        T_reservoir_limits_expanded_vector_degC = sort([T_init_vector_limits_degC T_reservoir_limits_degC]);
        [T_init_vector_degC_meshgrid,T_reservoir_degC_meshgrid] = meshgrid(T_init_vector_limits_degC,T_reservoir_limits_expanded_vector_degC);
        temporary_var=cat(2,T_init_vector_degC_meshgrid',T_reservoir_degC_meshgrid');
        all_combos_Tinit_Treservoir_C_acc = unique(reshape(temporary_var,[],2),'rows'); % Each row in this matrix contains all combinations of initial and reservoir temperatures
    elseif strcmp(temperature_permutations, 'extremes_only')
        middle_temperature = 25; % degC
        all_combos_Tinit_Treservoir_C_acc      = zeros(5,2); % Pre-allocate for the four extreme case & the "middle" case. Hence 5 rows
        all_combos_Tinit_Treservoir_C_acc(1,1) = min(T_init_vector_limits_degC);
        all_combos_Tinit_Treservoir_C_acc(2,1) = max(T_init_vector_limits_degC);
        all_combos_Tinit_Treservoir_C_acc(3,1) = middle_temperature;
        all_combos_Tinit_Treservoir_C_acc(4,1) = max(T_init_vector_limits_degC);
        all_combos_Tinit_Treservoir_C_acc(5,1) = min(T_init_vector_limits_degC);

        all_combos_Tinit_Treservoir_C_acc(1,2) = min(T_reservoir_limits_degC);
        all_combos_Tinit_Treservoir_C_acc(2,2) = min(T_reservoir_limits_degC);
        all_combos_Tinit_Treservoir_C_acc(3,2) = middle_temperature;
        all_combos_Tinit_Treservoir_C_acc(4,2) = max(T_reservoir_limits_degC);
        all_combos_Tinit_Treservoir_C_acc(5,2) = max(T_reservoir_limits_degC);
    end


    all_combos_Tinit_Treservoir_C_acc = all_combos_Tinit_Treservoir_C_acc(1:end-1,:);     % Remove unrealistc scenarios
    all_combos_Tinit_Treservoir_K_acc = T_abs_zero_K + all_combos_Tinit_Treservoir_C_acc; % Each row in this matrix contains all combinations of initial and reservoir temperatures
    no_of_temperature_combos          = size(all_combos_Tinit_Treservoir_K_acc,1);

    n_optimal_acc_given_Tcombo               = NaN(no_of_temperature_combos,1);
    P_wheels_cruise_vector                   = NaN(no_of_temperature_combos,1);
    exit_reason_flag_vector                  = cell(no_of_temperature_combos,1);
    optimal_nomin_cell_capacities_acc_vector = NaN(size(n_optimal_acc_given_Tcombo));
    sim_results_structs_acc                  = cell(no_of_temperature_combos,1);
    param_acc_structs                        = cell(no_of_temperature_combos,1);

    commandwindow;

    if strcmp(search_strat, 'binary')
        cut_point_layers = [];
    end

    %% Determine, for all temperature permutations, the optimal number of layers required for acceleration spec
    for idx_temperature_combo = 1:no_of_temperature_combos
        T_init_K      = all_combos_Tinit_Treservoir_K_acc(idx_temperature_combo,1);
        T_reservoir_K = all_combos_Tinit_Treservoir_K_acc(idx_temperature_combo,2);

        if idx_temperature_combo>1 && strcmp(search_strat, 'binary')
            if isnan(n_optimal_acc_given_Tcombo(idx_temperature_combo-1))
                cut_point_layers = [];
            else
                cut_point_layers = n_optimal_acc_given_Tcombo(idx_temperature_combo-1);
            end
        elseif idx_temperature_combo == 1 && strcmp(search_strat, 'linear')
            cut_point_layers = [];  % Required as a placeholder argument in the call to search_layers_acc_simplified
        end
        fprintf('\nNow simulating Temperature combo %2d of %2d, with T_init = %4.2f degC & T_reservoir = %4.2f degC\n\n',idx_temperature_combo,no_of_temperature_combos,T_init_K-T_abs_zero_K,T_reservoir_K-T_abs_zero_K);

        param_acc{1} = apply_params_for_layer_optim_fcn(SOC_init_acc_pct,OperatingMode,TemperatureEnabled,cutoffSOC,CutoverSOC,CutoffVoltage,CutoverVoltage,solidDiffusionMode,T_init_K,T_reservoir_K,T_max_K,sim_datalog_interval,axial_nodes,radial_nodes);
        [n_optimal_acc_given_Tcombo(idx_temperature_combo),sim_results_structs_acc{idx_temperature_combo},param_acc_structs{idx_temperature_combo},optimal_nomin_cell_capacities_acc_vector(idx_temperature_combo),P_wheels_cruise_vector(idx_temperature_combo), exit_reason_flag_vector{idx_temperature_combo}] ...
            = search_layers_acc(n_min,n_max,cut_point_layers,neg_to_pos_len_ratio,Mv_with_pack_overhead,Cd,Av,tf,v_base,eta_drivetrain,no_of_cells_in_pack,vf,param_acc, search_strat);

        if ~isnan(n_optimal_acc_given_Tcombo(idx_temperature_combo))    % Print the "completed" statement only if the result isn't a NaN - i.e. if a suitable layer count was determined
            fprintf('Completed optimal layer choice computation for acceleration run with T_init = %4.2f degC & T_reservoir = %4.2f degC\n\n',T_init_K-T_abs_zero_K, T_reservoir_K-T_abs_zero_K);
        end

        fprintf('\nFor all temperature combos up to & including this one (%2d of %2d), the optimal layer choice for acceleration run is % d\n',idx_temperature_combo,no_of_temperature_combos,max(n_optimal_acc_given_Tcombo(1:idx_temperature_combo)));
        fprintf('---------------------------------------------------------------------------------------------------------------------\n\n');
    end % End of temperature loop

end % End of compute_optim_layers_acc function

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: