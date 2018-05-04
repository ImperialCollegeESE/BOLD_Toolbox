% This is the main user-script for layer optimisation that calls all other functions and scripts as needed
% IMPORTANT: Please run this file from the root of the BOLD toolbox source codes folder (BOLD_src) as files/directories containing results will be written here.
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% Copyright year : 2018
% Last Updated   : Tue May  1 12:33:19 BST 2018
% License        : MIT License

% Key point 1:keeps the ratio L_neg/L_pos constant for all layer choices
% Key point 2:keeps the stack-thickness constant for all layer choices
clear; close all; clc; format long g; format compact;

%% User Setup: Required Before Every Run
platform                         = 'BEV';           % Use either 'BEV' or 'PHEV' (string data type) as argument to the vehicle specs function call
temperature_permutations_fastchg = 'all'; % Use either 'all' or 'extremes_only'. fastchg is presently set to "all" for thermal design space maps
temperature_resolution           = 5;               % degC, typically 5. Use 1 for high resolution thermal design space maps
% NOTE: the number of (T_init,T_sink) permutations is [number of elements in T_init_vector_limits_degC] times [number of unique elems in (T_init_vector_limits_degC union T_reservoir_limits_degC)]

search_strat = 'linear';   % Options are 'linear' or 'binary'. This choice is applied to both acceleration & fast charging cases
n_min        =  1;         % Numerically possible minimum number of layers is one (The user may set a more informed choice if desired)
n_max_choice = 'physical'; % Options are 'physical' or 'usable'. 'physical' represents the higher limit on n_max

% Charging power for the battery pack in Chevrolet Bolt BEV is 80kW, as per it's user's manual
% Furthermore, 80 kW is Level III fast charging as per: https://publish.illinois.edu/grainger-ceme/files/2014/06/CEME412_ChengGeorgiaTech1.pdf
P_batt_charging_kW = 50;  % kW

% workspace_Filename = 'Workspace50kWPHEV.mat';  % If user does not wish to save workspace (due to size), then please comment out this line and the last line in this script

%% Initial setup
set(0,'defaultaxesfontsize',12 ...
    ,'defaultaxeslinewidth',2 ...
    ,'defaultlinelinewidth',2.5 ...
    ,'defaultpatchlinewidth',2 ...
    ,'DefaultFigureWindowStyle','docked');

sim_datetime = datestr(datetime('now'));
sim_datetime = strrep(sim_datetime,':','_');
sim_datetime = strrep(sim_datetime,' ','_');

if exist('final_results_optimal_layer','dir')==0
    mkdir('final_results_optimal_layer');
end

%% Vehicle/Battery Specs & Physical num_layers
[Mv_with_pack_overhead,Cd,Av,acc_time_manuf,v_base,eta_drivetrain,no_of_cells_in_pack,...
    ~,vf_manuf, SOC_init_fast_chg_pct, CP_finish_pct_fast_chg] = vehicleSpecs(platform);  % Obtain vehicle specs (Mv = vehicle mass including pack overhead, Cd = coeff. of drag, Av = frontal area, acc_time = Acc. time specs, v_base = base speed [m/s], eta_drivetrain = drivetrain efficiency)
% Codes from this point onwards is based on LIONSIMBA toolbox. A modified version suitable for this layer optimisation task is provided here.
param{1}             = Parameters_init_suppliedSOC_pct(100); % LIONSIMBA dummy 'param' cell struct to extract the domain lengths and compute layer bounds
neg_to_pos_len_ratio = param{1}.len_n/param{1}.len_p;        % This is a critical constant that is a key idea of this paper
n_max_physical       = compute_max_layers_numerical_limit(neg_to_pos_len_ratio,param); % returns the max possible value of number of layers that can be used, keeping the same len/pos ratio, without domain thicknesses becoming negative.
clear param;

%% Layer optimisation for acceleration specs (0-60 mph)
acc_test_choice = 'worst_case'; % valid options are 'sae_derived' and 'worst_case'
fprintf('\nNow performing layer optimisation for acceleration specs ...\n\n');
T_init_vector_limits_acc_degC = 15:temperature_resolution:38; % Cell initial temperature range (as per SAE spec J1666 https://energy.gov/sites/prod/files/2015/04/f21/ntp002.pdf)
T_reservoir_limits_acc_degC   = 5:temperature_resolution:49;  % reservoir temperature range allowed during testing (reservoir range is much wider)
if temperature_resolution > 1                                 % Required to ensure the uppermost temperature values are included
    T_init_vector_limits_acc_degC = horzcat(T_init_vector_limits_acc_degC,38);
    T_reservoir_limits_acc_degC   = horzcat(T_reservoir_limits_acc_degC,49);
end

%% Determine acceleration params & lower layer bound
[vf, tf] = compute_tf_acc(acc_test_choice, acc_time_manuf, vf_manuf);
if strcmp(n_max_choice, 'physical')
    n_max = n_max_physical;
elseif strcmp(n_max_choice, 'usable')
    n_max = compute_nmax_usable(n_min, neg_to_pos_len_ratio, vf, tf, eta_drivetrain, Mv_with_pack_overhead, no_of_cells_in_pack, Cd, Av, v_base, n_max_physical);  % Compute the max. num layers, defined by the minimum of the ratio of power density to nominal capacity as num_layers increases
end

%% Determine optimum layer count for acceleration run
commandwindow;
[n_optimal_acc_given_Tcombo,sim_results_structs_acc,param_acc_structs, all_combos_Tinit_Tambient_K_acc, optimal_nomin_cell_capacities_acc_vector,P_wheels_cruise_vector, exit_reason_flag_vector_acc] = ...
    compute_optim_layers_acc(neg_to_pos_len_ratio,n_min,n_max,Mv_with_pack_overhead,Cd,Av,v_base,eta_drivetrain,no_of_cells_in_pack,...
    T_init_vector_limits_acc_degC, T_reservoir_limits_acc_degC, platform, vf, tf, search_strat);

%% Pretty-Printing acceleration run layer-optimisation results to the command window (in tabular form)
inital_temps_simulated_acc  = all_combos_Tinit_Tambient_K_acc(:,1) - 273.15;
ambient_temps_simulated_acc = all_combos_Tinit_Tambient_K_acc(:,2) - 273.15;
optimal_results_table_acc   = table(inital_temps_simulated_acc,ambient_temps_simulated_acc,ceil(n_optimal_acc_given_Tcombo),cellstr(num2str(optimal_nomin_cell_capacities_acc_vector,'%4.2f')), 'VariableNames',{'Init_Temp' 'Amb_Temp' 'Optimal_Layers_Acc','Nominal_Capacity_Ah'});
disp(sortrows(optimal_results_table_acc,3,'descend'));
writetable(optimal_results_table_acc,['final_results_optimal_layer/' 'acc_layer_results_' sim_datetime search_strat platform '.txt'],'Delimiter','\t');

[n_optimal_acc, idx_optimal_layer]   = max(n_optimal_acc_given_Tcombo);
param_acc_optimal                    = param_acc_structs{idx_optimal_layer}{1};
sim_results_acc_phase_optimal_layers = sim_results_structs_acc{idx_optimal_layer}{1};
cruise_distance                      = 1609.34; % [m] disctance to cruise at the final speed (after acceleration phase)
cruise_interval                      = cruise_distance/vf; % [sec] time-period in cruise phase
P_batt_cruise                        = -(P_wheels_cruise_vector(idx_optimal_layer)/eta_drivetrain); % Battery pack power (negative sign introduced here for adhering to LIONSIMBA conventions, i.e. discharge power is negative. Accelerating discharges the pack)
P_cell_cruise                        = P_batt_cruise/no_of_cells_in_pack; % Power to be supplied by a single cell
P_density_cruise                     = P_cell_cruise/param_acc_optimal{1}.overall_surface_area_for_given_layers;
sim_results_cruise_phase             = startSimulation(sim_results_acc_phase_optimal_layers.time{1}(end),sim_results_acc_phase_optimal_layers.time{1}(end)+cruise_interval,sim_results_acc_phase_optimal_layers.initialState,P_density_cruise,param_acc_optimal);
sim_results_acc_layer_optimisation   = [sim_results_acc_phase_optimal_layers;sim_results_cruise_phase];

%% Plot results of acceleration layer optimisation
sim_results_acc_combined_time_vector    = [sim_results_acc_layer_optimisation(1).time{1};sim_results_acc_layer_optimisation(2).time{1}];
sim_results_acc_combined_voltage_vector = [sim_results_acc_layer_optimisation(1).Voltage{1};sim_results_acc_layer_optimisation(2).Voltage{1}];
plot(sim_results_acc_combined_time_vector,sim_results_acc_combined_voltage_vector);grid on;title('Acceleration run with optimal layers');
ylabel('Terminal Voltage [V]'); xlabel('Time [sec]');

%% Determine optimum layer count for fast-charging
fprintf('\n---------------------------------------------------------------------------------------\n\n');
fprintf('\nNow performing layer optimisation for fast-charging ...\n\n');

T_init_vector_limits_fastchg_degC = 15:temperature_resolution:38; % Ref for init & reservoir temps for level III fast charge: https://energy.gov/sites/prod/files/2015/04/f21/ntp013.pdf
T_reservoir_limits_fastchg_degC   = 5:temperature_resolution:49;  % Reservoir temperature range allowed during charging. reservoir temps will be these AND those that are specified as init cell temps!
if temperature_resolution > 1                                     % Required to ensure the uppermost temperature values are included
    T_reservoir_limits_fastchg_degC = horzcat(T_reservoir_limits_fastchg_degC,49);
    T_init_vector_limits_fastchg_degC = horzcat(T_init_vector_limits_fastchg_degC,38);
end

[n_optimal_fastchg_given_Tcombo,sim_results_structs_fastchg,param_fastchg_structs, all_combos_Tinit_Tambient_K_fastchg, optimal_nomin_cell_capacities_vector_fastchg, exit_reason_flag_vector_fastchg,ix_saturation_vector,SOC_log_vector] = ...
    compute_optim_layers_fastchg(neg_to_pos_len_ratio,n_min,n_max,no_of_cells_in_pack,P_batt_charging_kW,T_init_vector_limits_fastchg_degC,T_reservoir_limits_fastchg_degC,...
    temperature_permutations_fastchg, platform, SOC_init_fast_chg_pct, CP_finish_pct_fast_chg, search_strat);

%% Pretty-Printing fast charge layer optimisation results to the command window (in tabular form)
inital_temps_simulated_fastchg  = all_combos_Tinit_Tambient_K_fastchg(:,1) - 273.15;
ambient_temps_simulated_fastchg = all_combos_Tinit_Tambient_K_fastchg(:,2) - 273.15;

usable_capacities_vector = cell(length(sim_results_structs_fastchg),1);
% figure(2);
% hold on;
legend_entries = num2str([inital_temps_simulated_fastchg ambient_temps_simulated_fastchg]);
for Tcombo = 1:length(sim_results_structs_fastchg)
    if iscell(sim_results_structs_fastchg{Tcombo})  % Only compute usable capacity if a successful layer was determined for this temperature combination
        surface_area = param_fastchg_structs{Tcombo}{1}{1}.overall_surface_area_for_given_layers; % Obtain surface area to convert curr_density to current in amps
        current      = sim_results_structs_fastchg{Tcombo}{1}.curr_density*surface_area;          % Coulomb-count current delivered on fast charge to the target SOC
        usable_capacities_vector{Tcombo} = trapz(sim_results_structs_fastchg{Tcombo}{1}.time{1}, (current/3600)); % Compute the Ah delivered for this layer configuration (i.e. USABLE or accessible capacity)
        % if strcmp(search_strat, 'linear')               % Only layers have been increased monotonically
        %    plot((n_min:n_max),SOC_log_vector{Tcombo}); % Visualise how the final SOC (upon saturation) at end of charging evolved with increasing number of layers
        % end
    else
        usable_capacities_vector{Tcombo} = NaN;         % If no layer was found, usable capacity is NaN
    end
end
if strcmp(search_strat, 'linear')
    legend(legend_entries);
end
hold off;

optimal_results_table = table(inital_temps_simulated_fastchg,ambient_temps_simulated_fastchg,ceil(n_optimal_fastchg_given_Tcombo),cellstr(num2str(optimal_nomin_cell_capacities_vector_fastchg,'%4.2f')), cellstr(num2str(cell2mat(usable_capacities_vector),'%4.2f')), 'VariableNames',{'Init_Temp' 'Amb_Temp' 'Optimal_Layers' 'Nominal_Capacity_Ah' 'Usable_Capacity_Ah'});
disp(sortrows(optimal_results_table,3,'descend'));
writetable(optimal_results_table,['final_results_optimal_layer/' 'fast_charge_layer_results_' sim_datetime search_strat platform num2str(P_batt_charging_kW) '.txt'],'Delimiter','\t');

%% Analyse results of charging
n_global_optimal_fastchg      = max(n_optimal_fastchg_given_Tcombo); % Optimum number of layers the greatest, since it covers the worst-case temperature scenario
indices_optimal_layer_fastchg = find(n_optimal_fastchg_given_Tcombo == n_global_optimal_fastchg); % Get indices of all temperature combinations that result in the max. num_layers

if ~isempty(indices_optimal_layer_fastchg)  % indices_optimal_layer_fastchg is empty if no optimal layer was determined, so only do the following if it's not empty..
    governing_temperature_combos_fastchg = all_combos_Tinit_Tambient_K_fastchg(indices_optimal_layer_fastchg,:);
    param_fastchg   = param_fastchg_structs(indices_optimal_layer_fastchg);       % Retrieve the parameter set correpsonding to the optimal layer
    param_fastchg   = param_fastchg{1};
    results_fastchg = sim_results_structs_fastchg(indices_optimal_layer_fastchg); % Retrieve the results correpsonding to the optimal layer
    results_fastchg = results_fastchg{1};
end

%% Save Workspace
% save(workspace_Filename, '-v7.3');

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t:
