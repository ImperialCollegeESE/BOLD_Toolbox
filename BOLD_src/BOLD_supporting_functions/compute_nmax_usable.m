function n_max_usable = compute_nmax_usable(n_min, neg_to_pos_len_ratio, vf, tf, eta_drivetrain, Mv_with_pack_overhead, no_of_cells_in_pack, Cd, Av, v_base, n_max_physical)
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

    param_100pct_soc = Parameters_init_suppliedSOC_pct(100);    % For computing nominal cell capacity.
    param_0pct_soc   = Parameters_init_suppliedSOC_pct(0);

    P_density_acc_log_local              = zeros(n_max_physical,1);
    P_density_acc_log_local(:,:)         = NaN;                     % Over-ride with NaNs since zeros have more possibility of messing with minimum determination
    nominal_cell_capacity_local_log      = zeros(n_max_physical,1);
    nominal_cell_capacity_local_log(:,:) = NaN;                     % Over-ride with NaNs since zeros have more possibility of messing with minimum determination

    num_layers_vector = n_min:n_max_physical;

    for num_layers = n_min:n_max_physical

        param{1} = Parameters_init_suppliedSOC_pct(100);   % Refresh/cleanse param struct for every num_layers

        % Re-compute parameters which are functions of number of layers
        [~, len_p, len_n] = compute_domain_thicknesses(num_layers, neg_to_pos_len_ratio, param);    % N.B. This param, & many others, must be replaced with param from apply_params_for_layer_optim_fcn if model calls are ever to be made
        param{1}.len_p    = len_p;
        param{1}.len_n    = len_n;

        param{1}.overall_surface_area_for_given_layers = num_layers*param{1}.surface_area_per_face_Northrop_cell;
        [param{1},~,~,~] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(num_layers, param{1}); % cell mass and Cp are appended to param struct

        %%%% Under acceleration %%%%
        Mv_local                 = Mv_with_pack_overhead + param{1}.mass_cell*no_of_cells_in_pack; % Entire vehicle mass with pack (neglecting rot. components)
        P_wheels_cruise_local    = compute_cruise_power(vf, Mv_local, Cd, Av);                     % Vehicle's power demand (watts) without considering power required to accelerate its mass.
        P_acc_mass_local         = 0.5*Mv_local*(v_base^2 + vf^2)/tf;                              % Power required solely for accelerating the vehicle's mass
        P_wheels_acc_total_local = P_wheels_cruise_local + P_acc_mass_local;                       % Total power demand during acceleration phase (0-60mph)

        P_batt_acc_local    = -(P_wheels_acc_total_local/eta_drivetrain); % Battery pack power (negative sign introduced here for adhering to LIONSIMBA conventions, i.e. discharge power is negative. Accelerating discharges the pack)
        P_cell_acc_local    = P_batt_acc_local/no_of_cells_in_pack;       % Power to be supplied by a single cell
        P_density_acc_local = P_cell_acc_local/param{1}.overall_surface_area_for_given_layers;
        P_density_acc_log_local(num_layers) = P_density_acc_local;

        %%%% Nominal capacity %%%%
        nominal_cell_capacity_local = compute_capacity_for_layer_fcn(num_layers,param_100pct_soc,param_0pct_soc,len_n);
        nominal_cell_capacity_local_log(num_layers) = nominal_cell_capacity_local;

    end     % End of while loop

    %%%% Ratio %%%%
    Pdensity_NomCap_ratio           = abs(P_density_acc_log_local./nominal_cell_capacity_local_log);                                                                                                                            % abs since acceleration power is negative
    min_Pdensity_NomCap_ratio_index = find(Pdensity_NomCap_ratio == min(Pdensity_NomCap_ratio(:)));                                                                                                                 =  = min(Pdensity_NomCap_ratio(:)));
    index_for_num_layers_vec        = (min_Pdensity_NomCap_ratio_index - n_min + 1);                                                                                                                                            % Index needs to be negatively offset by the n_min value + 1
    n_max_usable                    = num_layers_vector(index_for_num_layers_vec);                                                                                                                                              % Number of layers at which binary search algorithm ceases to be valid, since ratio changes sign here

end     % End of function

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: