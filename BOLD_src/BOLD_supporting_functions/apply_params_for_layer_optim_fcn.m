function param = apply_params_for_layer_optim_fcn(init_cell_soc_percent,OperatingMode,TemperatureEnabled,cutoffSOC,CutoverSOC,CutoffVoltage,CutoverVoltage,solidDiffusionMode,T_init,T_ambient,Tmax,sim_datalog_interval,axial_nodes, radial_nodes)
% Import LIONIMBA params and over-ride relevant quantities (other than the thermal model, which will be done within the main script)
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

    param                 = Parameters_init_suppliedSOC_pct(init_cell_soc_percent); % Load up the parameters from the Parameters_init_suppliedSOC_pct file
    param.PrintHeaderInfo = 0; % Suppress printing of headers for speed
    param.Scope           = 0; % turn off the scope for speed

    param.OperatingMode      = OperatingMode;
    param.TemperatureEnabled = TemperatureEnabled;

    param.CutoffSOC      = cutoffSOC;       % [%] Cut-off SOC percentage
    param.CutoverSOC     = CutoverSOC;      % [%] Cut-over SOC percentage
    param.CutoffVoltage  = CutoffVoltage;   % [V] Define lower cutoff voltage
    param.CutoverVoltage = CutoverVoltage;  % [V] Define upper cutover voltage

    param.SolidPhaseDiffusion = solidDiffusionMode;

    param.T_init = T_init;    % K
    param.Tref   = T_ambient; % K
    param.Tmax   = Tmax;

    param.sim_datalog_interval = sim_datalog_interval;

    %% Refine node densities
    % param.Nal  = 10; % Not relevant if using a lumped thermal model % Number of control volumes used for discretising the aluminium current collector domain in the axial (through-thickness) direction
    param.Np   = axial_nodes; % Number of control volumes used for discretising the positive electrode domain in the axial (through-thickness) direction
    param.Ns   = axial_nodes; % Number of control volumes used for discretising the separator domain in the axial (through-thickness) direction
    param.Nn   = axial_nodes; % Number of control volumes used for discretising the negative electrode domain in the axial (through-thickness) direction
    % param.Ncu  = 20; % Not relevant if using a lumped thermal model % Number of control volumes used for discretising the copper current collector domain in the axial (through-thickness) direction

    param.Nr_p = radial_nodes; % Number of control volume (shells) used for discretising each cathode particle in the radial/spherical direction (P2D direction)
    param.Nr_n = radial_nodes; % Number of control volume (shells) used for discretising each anode particle in the radial/spherical directi

end
% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: