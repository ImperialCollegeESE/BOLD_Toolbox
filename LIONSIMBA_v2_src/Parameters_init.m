function param= Parameters_init
% PARAMETERS_INIT Defines the parameters used in simulation. To change the cell/simulator parameters, modify the values here.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to marcello.torchio01@ateneopv.it

% Copyright 2015:   Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%                   Bhushan Gopaluni, University of British Columbia
%                   Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init_cell_soc_percent = 85.51; % Cell-SOC [%] at start of simulation for compatibility with existing LIONSIMBA code (<=1.023)

param.OperatingMode = 1;

% Valid options are:
%                       1 - Constant input current density(default mode) between the initial time
%                           and final time.
%                       2 - Constant input power density between the initial time
%                           and final time.
%                       3 - Potentiostatic charge. In this case the current
%                           is considered a variable and the charge is
%                           carried out at constant potential.
%                       4 - Variable current density profile described in the
%                           getInputCurrent script as a function of t
%                       5 - Variable power density profile described in the
%                           getPowerCurrent script as a function of t

% choice of (lumped) thermal model
param.lumped_thermal_version = 1;

% Valid options are:
%                       1 - Lumped thermal model with heat generation from
%                           eta*current only
%                       2 - Lumped thermal model with heat generation from
%                           both entropy & eta*current

% Option to suppress command window output from 'startSimulation'
param.suppress_status_prints = 0;   %  1 = Suppressed output. 0 = output as normal

% Stoichiometry Limits
param.theta_max_pos = 0.49550;  % at 100% cell SOC
param.theta_max_neg = 0.85510;  % at 100% cell SOC
param.theta_min_pos = 0.99174;  % at 0% cell SOC
param.theta_min_neg = 0.01429;  % at 0% cell SOC

%% Universal Constants
% Faraday Constant  [C/mol]
param.F         = 96487;
% Gas constant      [J / (mol K)]
param.R         = 8.314;

%% Cell Exterior Geometry (Pouch Cell)
param.pouch_length = 332.74e-3;         % [m] (long dimension)  Length of cell pouch
param.pouch_width  = 99.06e-3;          % [m] (short dimension) Width of cell pouch
tab_width = 40e-3;       % presently used for current collection as well as cooling (surface cooling will result in smaller tabs)
tab_length = 0.75*param.pouch_width;     % about 70%-80% of the pouch's dimension, where the length (L) represents shorter pouch dimension
param.tab_area = 2*tab_width*tab_length;  % there are 2 tabs available for convecting heat out (used only in tab-cooling approach)

%% Sections thickness    [m]
% Aluminium current collector
param.len_al= 15e-6;
% Positive Electrode
param.len_p = 72e-6;
% Separator
param.len_s = 25e-6;
% Negative Electrode
param.len_n = 88e-6;
% Copper current collector
param.len_cu= 10e-6;
% Pouch Wrapper
param.len_pouch = 160e-6; % ref: "Li-Ion Pouch Cells for Vehicle Applications—Studies of Water Transmission and Packing Materials", Pontus Svens, Maria Hellqvist Kjell, Carl Tengstedt, Göran Flodberg and Göran Lindbergh, Energies 2013, 6, 400-410; doi:10.3390/en6010400

%% Cell Electrical/Internal Geometry/Design Aspects (extrapolated from simulations from LIONSIMBA<=1.023 and compared to standard literature)
% param.i_1C_density_Northrop_cell below is the discharge current density required in order to deplete the (Northrop) cell beginning
% at 100% SOC down to an (arbitrary) cut-off of 2.7V for the LCO cell with parameters described here (Northrop cell).
% This represents an indirect measure of capacity of one electrochemical layer (Al-Pos-Neg-Sep-Cu combo)(that is simulated by the original Newman model as published in typical literature.)

param.i_1C_density = 29.23; % [A/m^2]
param.no_of_layers_Northrop_cell = 49;     % Assumed no. of layers for LIONSIMBA's reference cell (Northrop cell), obtained by hand-computing the number of layers that can be fit inside a 10mm pouch.
param.t_stack = param.no_of_layers_Northrop_cell*(param.len_p + param.len_s + param.len_n) + (ceil(0.5*(param.no_of_layers_Northrop_cell + 1))*param.len_cu) + ceil(0.5*param.no_of_layers_Northrop_cell)*param.len_al; % length/thickness inside pouch available for filling up with unit cells (i.e. the stack of layers)
assumed_cell_capacity_Ah = 60; % Assumed Ah of a cell constructed with multiple layers (each layer having the parameters given in this file).

% NOTE: The code does not simulate a cell of the capacity given above. In fact, LIONSIMBA only simulates just one electrochemical layer.
% The variable above is used only for calculating the overall surface area of all electrochemical layers (used in lumped thermal params and potentially useful in user-level scripts)
I_1C_cell_amps = assumed_cell_capacity_Ah; % By definition of C-rate. (i.e. a 60Ah cell has a 1C current of 60A)

% The variable computed below represents the area in the plane that is perpendicular to the through-thickness direction in a typical 1D discretisation of a standard Newman model)
param.overall_surface_area_for_given_layers = I_1C_cell_amps/param.i_1C_density;  % [m^2] overall surface area provided by all layers

%% Thermal conductivities [ W / (m K) ]

% Aluminium current collector
param.Lambda_al = 237;  % From material datasheets and standard scientific tables
% Positive electrode
param.Lambda_p  = 1.85; % Ref: Table 1 from Thermal modeling of cylindrical lithium ion battery during discharge cycle", % Dong Hyup Jeon, Seung Man Baek, Energy Conversion and Management, 2011
% Separator
param.Lambda_s  = 0.35; % around 0.1 to 1. A typical value of 0.35 is assumed, since we do not know the separator material itself
% Negative Electrode
param.Lambda_n  = 5;    % Ref: Table 1 from Thermal modeling of cylindrical lithium ion battery during discharge cycle", % Dong Hyup Jeon, Seung Man Baek, Energy Conversion and Management, 2011
% Copper current collector
param.Lambda_cu = 401;  % From material datasheets and standard scientific tables

%% Electrolyte diffusion coefficients [m^2 / s]

% Positive domain
param.Dp = 3.2227e-10;  % USED ONLY FOR NON THERMALLY COUPLED SIMULATION. ASSUMES 25 DEGREES C
% Separator
param.Ds = 3.2227e-10;  % USED ONLY FOR NON THERMALLY COUPLED SIMULATION. ASSUMES 25 DEGREES C
% Negative domain
param.Dn = 3.2227e-10;  % USED ONLY FOR NON THERMALLY COUPLED SIMULATION. ASSUMES 25 DEGREES C

%% Density [kg / m^3 ]
% Aluminium current collector
param.rho_al = 2700; % At room temp. From material datasheets and standard scientific tables
% Positive electrode
param.rho_p  = 2291.62; % Ref: Table 1 from Thermal modeling of cylindrical lithium ion battery during discharge cycle", Dong Hyup Jeon, Seung Man Baek, Energy Conversion and Management, 2011
% Separator
param.rho_s  = 1100; % We do not know exactly the separator material used. Common separator materials have a density of 1070 - 1200 kg/m^3  (1070, 1100, 1200 etc.). Using a median value
% Negative electrode
param.rho_n  = 2660; % Table 1 from 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior Using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173 (depends on the exact carbonaceous material used)
% Copper current collector
param.rho_cu = 8960; % At room temp. From material datasheets and standard scientific tables
% LiPF6 electrolyte
param.rho_LiPF6 = 1290; % Table 1 from 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior Using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173, also supported by 'Thermal analysis of a cylindrical lithium-ion battery', Xiongwen Zhang, Electrochimica Acta,  56 (2011) 1246–1255, Table 3
% Pouch Material
param.rho_pouch = 1150; % 'Modeling for the scale-up of a lithium-ion polymer battery',Ui Seong Kim, Chee Burm Shin, Chi-Su Kim, Journal of Power Sources, 2008
% Filler/Binder
param.rho_pvdf = 1750;  % Table 1 from 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173

%% Temperature Settings
param.Tref   = 25 + 273.15; % Environment (ambient) temperature [K]

% ref: 'Reciprocating air flow for Li-ion battery thermal management to improve temperature uniformity" , Rajib Mahamud, Chanwoo Park, Journal of Power Sources, 196 (2011) 5685–5696
param.Tmax   = 55 + 273.15; % Absolute maximum permissible temperature [K] Upper limit on temperature, at ANY point in ANY cell (in the pack) during operation

%% Specific heat capacities [ J / (kg K) ]
param.Cpal    =   903.0;  % Aluminium current collector
param.Cpp     =  1269.2;  % Positive Electrode
param.Cps     =  1978.2;  % Separator
param.Cpn     =  1437.4;  % Negative Electrode
param.Cpcu    =   385.0;  % Copper current collector
param.CpLiPF6 =  2055.1;  % Electrolyte (from S. C. Chen, C. C. Wan, Y. Y. Wang, Thermal analysis of lithium-ion batteries, Journal of Power Sources 140 (1) (2005) 111124. doi: DOI: 10.1016/j.jpowsour.2004.05.064. URL http://www.sciencedirect.com/science/article/B6TH1- 4DF49BM-2/2/a6c383e4ffa32dba3641a1e863344ae9)
% param.CpLiPF6 =  134.1; % Electrolyte (too low, invalid)
% Assumption: Ignoring Cp of binder/filler since they are negligible in content
% furthermore, the exterior pouch is also ignored in Cp calculations (but it is accounted for in mass calculations)
param.Cppouch = 1464.8; % Weighted calculation based on the constituents of the pouch material

%% Current collector conductivities [S/m]
param.sig_al = 3.55e7;
param.sig_cu = 5.96e7;

%% Electrolyte Porosity indices (i.e. electrolyte volume fractions)
% Positive electrode
param.eps_p     = 0.385;
% Separator
param.eps_s     = 0.724;
% Negative electrode
param.eps_n     = 0.485;


%% Volume fraction
param.eps_fi    = [0.025;0;0.0326];

%% Bruggeman coefficients

% Positive electrode
param.brugg_p   = 4;
% Separator
param.brugg_s   = 4;
% Negative electrode
param.brugg_n   = 4;

%% Solid diffusion coefficients [m^2 / s]

% Positive electrode
param.Dps       = 1e-14;
% Negative electrode
param.Dns       = 3.9e-14;

%% Particle surface area [m^2 / m^3]
% Positive electrode
a_p       = 885000;
% Separator
a_s       = 0;
%Negative electrode
a_n       = 723600;
% Do not remove. It is used in the code
param.a_i       = [a_p;a_s;a_n];

%% Transference number - Not available for positive/negative electrode
param.tplus     = 0.364;

%% Reaction rate constants [ m ^ 2.5 / (mol^0.5 s ) ]
% Positive electrode
param.k_p       = 2.334e-11;
% Separator
param.k_s       = 0;
% Negative electrode
param.k_n       = 5.031e-11;

%% Heat exchange coefficient [W/m^2 K]
param.hcell     = 1;                    % Used only for 1D thermal model in through-thickness axis

%% Maximum concentration of Li-ions in the solid phase [ mol/m^3 ]
% Positive electrode
param.cs_maxp   = 51554;
% Separator
param.cs_maxs   = 0;
% Negative electrode
param.cs_maxn   = 30555;

% Solid particle radius [m] - It's equal for the positive and negative electrode
param.Rp_p     = 2e-6;
param.Rp_n     = 2e-6;

% Solid phase conductivities (S/m)
param.sig    = [100;... % Positive electrode
    0;...   % Separator
    100     % Negative electrode
    ];

param.vol_fraction_solidphase = (1 - [param.eps_p;param.eps_s;param.eps_n] - param.eps_fi);
param.vol_fraction_solidphase(2) = 0;   % No solid phase material in separator
% Effective solid phase conductivities (S/m)
param.sig_eff = param.sig.*param.vol_fraction_solidphase;

%% Lumped Thermal Parameters
[param,surface_area_per_face_for_49_layers,~,~] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(param.no_of_layers_Northrop_cell,param); % computes the cell's mass & Cp_avg and appends them to param struct (to be used in lumped thermal model calcs)
param.surface_area_per_face_Northrop_cell = surface_area_per_face_for_49_layers;   % This will be the surface area per face specific to the 49 layer Northrop cell mapped to a 10mm thick pouch of spefic WXL
param.h_lumped = 150;    % [W/(m^2 K)] IDA seems to be very sensitive to this parameter during charging. Please over-write this in user-scripts with a suitable value (this value is highly dependent upon the specific cooling used)
clear surface_area_per_face_for_49_layers;

%% Activation Energy for Temperature Dependent Solid Phase Diffusion [ J / mol ]

% Positive electrode
param.EaDps  = 5000;

% Negative Electrode
param.EaDns  = 5000;

%% Activation Energy for Temperature Dependent Reaction Constant [ J / mol ]

% Positive electrode
param.Eakip  = 5000;

% Negative electrode
param.Eakin  = 5000;


%% Initial conditions

% Electrolyte Li-ions initial concentration [mol/m^3]
param.ce_init = 1000;

% Initial temperature of the cell [K]
param.T_init = 298.15;


%% Simulator parameters

% Select the model used for approximate the solid phase diffusion
% Allowed values are:
%                       1 - Parabolic approximation (two parameters model)
%
%                       2 - Higher-order polynomial (three parameters
%                       model)
%
%                       3 - Full order model
%

param.SolidPhaseDiffusion = 1;

% Enable or disable temperature dynamics

param.TemperatureEnabled = 1;

% Integration step [s]
param.sim_datalog_interval = 0.5; % interval for logging data (although IDA is an adaptive solver, at the end of simulation, the user may like output variables at a specific time interval)

% Cutoff voltage [V]
param.CutoffVoltage = 2.5;

% Cutover voltage [V]
param.CutoverVoltage = 4.3;

% Cutoff SOC [%]
param.CutoffSOC = 0;

% Cutover SOC [%]
param.CutoverSOC = 100;

% Number of control volumes used for discretising the aluminium current collector domain in the axial (through-thickness) direction
param.Nal   = 10; % Not relevant if using a lumped thermal model

% Number of control volumes used for discretising the positive electrode domain in the axial (through-thickness) direction
param.Np    = 10;

% Number of control volumes used for discretising the separator domain in the axial (through-thickness) direction
param.Ns    = 10;

% Number of control volumes used for discretising the negative electrode domain in the axial (through-thickness) direction
param.Nn    = 10;

% Number of control volumes used for discretising the copper current collector domain in the axial (through-thickness) direction
param.Ncu   = 10; % Not relevant if using a lumped thermal model


% If the full diffusion model (Fick's law) is selected, the below two parameters define the number of discretization points inside the solid particles.
% Number of control volume (shells) used for discretising each cathode particle in the radial/spherical direction (P2D direction)
param.Nr_p = 10;

% Number of control volume (shells) used for discretising each anode particle in the radial/spherical directi
param.Nr_n = 10;

% Initial concentration of Li-ions in the solid phase [mol/m^3]
param.init_cell_soc = init_cell_soc_percent/100; % convert to a fraction between 0 and 1
% Positive electrode
param.cs_p_init = ((param.init_cell_soc*(param.theta_max_pos-param.theta_min_pos) + param.theta_min_pos))*param.cs_maxp;

% Negative electrode
param.cs_n_init = ((param.init_cell_soc*(param.theta_max_neg-param.theta_min_neg) + param.theta_min_neg))*param.cs_maxn;
param.cs_neg_saturation = (0.01*param.CutoverSOC*(param.theta_max_neg-param.theta_min_neg) + param.theta_min_neg)*param.cs_maxn;
param.enable_csneg_Saturation_limit = 0; % This parameter, when set to 1, enforces simulation termination when surface concentration of any node in the neg electrode hits the saturation value (set in parameter above)
param.cs_sat_thresh = 1.0;               % Threshold cs fraction for fast-charging algorithm

% Enable or disable the scope in the matlab command line
param.Scope = 1;

% Enable or disable the printing of header information
param.PrintHeaderInfo = 1;

%% External functions

% This field can be used as an extra structure and it is passed to all the
% external scripts.
param.extraData = [];

% Define the name of the external function that has to be called to compute
% the value of applied current. This function is called during the integration process.

param.CurrentDensityFunction = @getInputCurrentDensity; % contains the handle to the variable current density profile (external function file)
param.PowerDensityFunction = @getInputPowerDensity;     % contains the handle to the variable power density profile (external function file)

% Define the name of the external function used to compute physical and
% transport properties of the materials during simulation. Please refer to
% the existing functions to get insight for custom implementations.

% Electrolyte diffusion coefficients
param.ElectrolyteDiffusionFunction          = @electrolyteDiffusionCoefficients;
% Electrolyte conductivity coefficients
param.ElectrolyteConductivityFunction       = @electrolyteConductivity;
% Open circuit potential
param.OpenCircuitPotentialFunction          = @openCircuitPotential;
% Solid phase diffusion coefficient
param.SolidDiffusionCoefficientsFunction    = @solidPhaseDiffusionCoefficients;
% Reaction rates
param.ReactionRatesFunction                 = @reactionRates;

% If a function handle is provided, the function is called after each
% integration step. All the states relative to the current integration step
% are provided besides extraData structure and timing information. See
% socEstimator.m as example
param.SOC_estimation_function = @socEstimator;

%% Potentiostatic mode
% This value (applicable only if the OperatingMode flag is set to 3) is used to control the battery in a potentiostatic manner.
param.V_reference = 4.2;

%% Tolerances
% Integrator (IDA) tolerances
param.AbsTol = 1e-6;
param.RelTol = 1e-6;

%% Ageing parameters (TESTING PURPOSES, BETA VERSION)

param.EnableAgeing = 0;

% Initial SEI resistance value [Ohm m^2]
param.R_SEI     = 0.01;
%Molar weight                               [kg/mol]
%ATTENTION: In Development of First Principles Capacity Fade Model
%for Li-Ion Cells, Ramadass et al. the measurement unit of M_p is wrong as
%well as the number itself. Please refer to Review of models for predicting
%the cycling performance of lithium ion batterise, Santhanagopalan et al.
param.M_n               = 73e-3;
% Admittance                                [S/m]
param.k_n_aging         = 3.79e-7;
% Side reaction current density             [A/m^2]
param.i_0_jside         = 0.80e-10;
% Open circuit voltage for side reaction    [V]
param.Uref_s            = 0.4;
% 1C current for the particular chemistry [A/m^2]
param.I1C               = 29.23;
% Weight used in the aging dynamics
param.w		= 2;
% for lumped thermal model dynamics, a cell weight calculation is given
% below. Perhaps the user might want to use that.
% param.w = param.mass_cell; % the cell mass calculation performed in this section might be useful to replace the default value of 2 used in the ageing dynamics section above

% Set to 1 if the user wants to use tha Jacobian matrix during
% calculations.
param.UseJacobian       = 1;

% This value, if set, represents the Jacobian function used from the
% integrator. It has to be a class 'Function' of the CasADi package. If
% provided, with UseJacobian=1, it will be used for speed up the
% integration process. If not provided, with UseJacobian=1, the code will
% compute the Jacobian on its own.
param.JacobianFunction = [];


% Type of the DAE system returned by LIONSIMBA batteryModel.m script. This function is under development.
% Admitted values are:
%                       1 - The equations are returned in an analytical
%                       form, written as implicit DAEs, i.e.
%
%                       x_dot - f(x,z) | Time differential equations
%                       z - g(x,z)     | Algebraic equations
%
%                       2 - The equations are returned in an analytical
%                       form, where time differential equations are written
%                       in an explicit form, i.e.
%
%                       f(x,z)         | Time differential equations
%                       z - g(x,z)     | Algebraic equations
%
param.daeFormulation = 1;

end
