function [Mv_with_pack_overhead,Cd,Av,acc_time_manuf,v_base,eta_drivetrain,no_of_cells_in_pack, regen_braking_fraction, vf_manuf, SOC_init_fast_chg_pct, CP_finish_pct_fast_chg] = vehicleSpecs(platform)
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

    pack_mass_chev_bolt    = 435.0;                                        % kg, (not directly used in layer opt calcs) source: Wikipedia and URL: http://media.chevrolet.com/media/us/en/chevrolet/vehicles/bolt-ev/2017.tab1.html
    no_of_pax              = 2;                                            % Number of passengers for acceleration testing
    mass_per_pax           = 75.3;                                         % kg (as per ETA-NTP002 standard, Revision 3, Effective December 1, 2004)
    cargo_wt_in_trunk      = 0;                                            % kg (weight of cargo, as per above standard)
    payload                = (no_of_pax*mass_per_pax) + cargo_wt_in_trunk; % This is the total payload to be carried in the vehicle
    Cd                     = 0.308;                                        % source: http://www.hybridcars.com/2017-chevy-bolt-ev-is-less-of-a-drag-than-originally-believed/
    vf_manuf               = 60;                                           % [mph] manufacturer spec for acceleration
    vehicle_height_inches  = 62.9;                                         % http://www.chevrolet.com/bolt-ev-electric-vehicle/specs/trims.html
    vehicle_height         = 0.0254*vehicle_height_inches;                 % inches to metres
    track_width_inches     = 59.1;                                         % source:Chevrolet Bolt EV US Website
    vehicle_track_width    = 0.0254*track_width_inches;                    % inches to metres
    Av                     = vehicle_height*vehicle_track_width;           % Very conservative estimation, without considering specific body shape/streamlining etc.
    acc_time_manuf         = 6.5;                                          % 0 - 60mph time, seconds for Chevrolet Bolt: http://www.chevrolet.com/bolt-ev-electric-vehicle-2
    v_base_mph             = 30.0;                                         % mph (base speed)
    v_base                 = 0.44704*v_base_mph;                           % metres/sec
    regen_braking_fraction = 0.85;                                         % Fraction of total braking power system can accept
    eta_drivetrain         = 0.75;                                         % overall drivetrain efficiency (simplified)
    cells_in_series        = 96;                                           % no. of series-connected cells in pack

    param{1}                   = Parameters_init_suppliedSOC_pct(100);     % 'dummy' param for lumped mass calculation for baseline cell in LIONSIMBA toolbox
    [~,~,mass_Northrop_cell,~] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(param{1}.no_of_layers_Northrop_cell,param{1});

    cells_in_parallel_BEV            = 3;  % no. of parallel-connected cells in pack
    no_of_cells_in_pack_BEV          = cells_in_series*cells_in_parallel_BEV;  % total number of cells in pack
    pack_mass_chev_bolt_overhead_BEV = pack_mass_chev_bolt - mass_Northrop_cell*no_of_cells_in_pack_BEV;  % Computing pack mass less the mass of the cells, for the default Bolt configuraiton (kg)
    CP_finish_pct_fast_chg           = 80; % Percent SOC that the cell should be at by end of CP phase, prior to pulsing

    if strcmp(platform,'BEV')
        curb_mass                 = 1625; % kg ,source: Chevrolet Bolt EV US Website
        Mv_with_pack_overhead     = curb_mass + payload + pack_mass_chev_bolt_overhead_BEV - pack_mass_chev_bolt; % Overall vehicle mass (neglecting components with rot. inertia) for all calculations
        no_of_cells_in_pack       = no_of_cells_in_pack_BEV;  % For function return
        SOC_init_fast_chg_pct_BEV = 20;   % Percent, for BEV, in accordance with platform-specific SOC limits
        SOC_init_fast_chg_pct     = SOC_init_fast_chg_pct_BEV;
    elseif strcmp(platform,'PHEV')
        curb_mass                         = (1625+98); % kg, 98 kg Ecotec engine mass added for PHEV
        cells_in_parallel_PHEV            = 1;         % no. of parallel-connected cells in pack
        no_of_cells_in_pack_PHEV          = cells_in_series*cells_in_parallel_PHEV;  % total number of cells in pack
        no_of_cells_in_pack               = no_of_cells_in_pack_PHEV; % For function return
        pack_mass_chev_bolt_overhead_PHEV = pack_mass_chev_bolt_overhead_BEV*(cells_in_parallel_PHEV/cells_in_parallel_BEV);  % Computing pack mass less the mass of the cells (kg)
        Mv_with_pack_overhead             = curb_mass + payload + pack_mass_chev_bolt_overhead_PHEV - pack_mass_chev_bolt; % Overall vehicle mass (neglecting components with rot. inertia) for all calculations
        SOC_init_fast_chg_pct_PHEV        = 30; % Percent, for PHEV, in accordance with platform-specific SOC limits
        SOC_init_fast_chg_pct             = SOC_init_fast_chg_pct_PHEV;
    else
        error('\nVehicle type must be either BEV or PHEV !....\n');
    end

end

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: