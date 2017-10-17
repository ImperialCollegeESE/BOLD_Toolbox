function results = startSimulation(t0,tf,initialState,input_density,startParameters)
% STARTSIMULATION  Starts the simulation of the Li-ion battery.
%
%   RESULTS = STARTSIMULATION(t0,tf,initialStates,INPUT_DENSITY,startParameters)
%   Starts the simulation of the Li-ion Cell with LIONSIMBA.
%
%   Input:
%       - t0 : initial integration time
%       - tf : final integration time
%       - initialStates : structure containing data for initializing the
%       states of the battery.
%       - INPUT_DENSITY  : Applied current/power density. If negative, the battery gets
%       discharged. If positive, the battery gets charged.
%       - startParameters : if provided, it defines the cell array containing the parameters
%       structures to be used in the simulation. Every single structure has to be obtained through the Parameters_init script.
%		If a cell array of N parameters structures is used, the simulator will perform a simulation with N cells in series.
%		If a cell of 1 parameters structure is used, then a single cell will be simulated.
%
%   Output:
%       - results : structure containing the solution of the dependent
%       variables. The stored results have as many rows as time instants
%       and columns as number of discrete volumes. If multiple cells are
%       simulated, the index i is used to access to the data of the i-th
%       cell.
%
%     results.Phis{i}:                      Solid Phase potential
%     results.Phie{i}:                      Liquid Phase potential
%     results.ce{i}:                        Electrolyte concentration
%     results.cs_surface{i}:                Electrode surface concentration
%     results.cs_average{i}:                Electrode average concentration
%     results.time{i}:                      Interpolated simulation time
%     results.int_internal_time{i}:         Integrator time steps
%     results.ionic_flux{i}:                Ionic flux
%     results.side_reaction_flux{i}:        Side reaction flux
%     results.SOC{i}:                       State of Charge
%     results.SOC_estimated{i}:             State of Charge estimate according to the
%                                           user-defined function
%     results.Voltage{i}:                   Cell voltage
%     results.Temperature{i}:               Cell Temperature
%     results.Qrev{i}:                      Reversible heat generation rate
%     results.Qrxn{i}:                      Reaction heat generation rate
%     results.Qohm{i}:                      Ohmic heat generation rate
%     results.film{i}:                      Side reaction film resistance
%     results.R_int{i}:                     Internal resistance
%     results.Up{i}:                        Cathode open circuit potential
%     results.Un{i}:                        Anode open circuit potential
%     results.etap{i}:                      Cathode overpotential
%     results.etan{i}:                      Anode overpotential
%     results.parameters{i}:                Parameters used for the simulation
%     results.JacobianFun:                  If evaluated, contains the
%                                           Jacobian matrix

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%                   Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% code infrastructure required to facilitate parallel simulations in an outer parfor loop with different parameter sets
% if iscell(startParameters)
%     temporary_variable = startParameters{1};
%     clear startParameters;
%     startParameters{1}=temporary_variable;
% else
%     temporary_variable = cell(1,1);
%     temporary_variable{1} = startParameters;
%     clear startParameters;
%     startParameters=temporary_variable;
% end

% Version of LIONSIMBA
version       = '2.0';

if(isempty(startParameters))
    % Load battery's parameters if not provided by the user
    param{1} = Parameters_init;
else
    % Use user provided parameters
    param = startParameters;
end

% Validate the input current density/power_density values
if(param{1}.OperatingMode~=3)
    if(~isreal(input_density) || isnan(input_density) || isinf(input_density) || isempty(input_density))
        error('The input current/power densities provided by user is complex valued or is NaN. Please check the values and restart simulation.')
    end
end

% Check for environmental tool availability
checkEnvironment(param,nargin);

% If enabled, print the header information
if(param{1}.PrintHeaderInfo==1)
    headerInfo(version)
end

% If everything is ok, let's start to simulate.
results = mainCore(t0,tf,initialState,input_density,param);
end

function results = mainCore(t0,tf,initialState,input_density,param)

% Store the original parameters structure in order to return it at the
% end of simulations.
param_original  = param;

% Get the total number of cells that have to be simulated.
n_cells         = length(param);

% Check if more cells are simulated when potentiostatic conditions are required.
if(param{1}.OperatingMode==3 && n_cells~=1)
    clc;
    error('!!!ERROR!!! -- Potentiostatic simulations are only possible with single cell packs -- !!!ERROR!!!')
end

% Check if the initial state structure is given
[Y0_existence,YP0_existence,Y0,YP0] = checkInitialStates(initialState);

% Switch among the selected operating modes for defining a suitable getCurrentDensity/getPowerDensity
% function. If multiple cells are being simulated, the variable current/power profile
% or its constant value are retreived from the first element of the
% parameters structures. This is valid because in a series string, all cells carry the same current.
switch(param{1}.OperatingMode)
    case 1
        param{1}.getCurrentDensity = @(t,t0,tf,extra)input_density;
    case 2
        param{1}.getPowerDensity = @(t,t0,tf,extra)input_density;
    case 3
        param{1}.getCurrentDensity = @(t,t0,tf,extra)0; % dummy values, will be over-written by solver
        param{1}.getPowerDensity = @(t,t0,tf,extra)0;   % unsure of this?
    case 4
        param{1}.getCurrentDensity = param{1}.CurrentDensityFunction; % external, i.e. user-supplied function values from external files loaded in parameterisation file
    case 5
        param{1}.getPowerDensity = param{1}.PowerDensityFunction; % external, i.e. user-supplied function values from external files loaded in parameterisation file
    otherwise
        error('Operating mode not supported');
end

% Define absolute and relative tolerances. If more cells are required,
% these values are taken from the first parameters structure.
opt.AbsTol      = param{1}.AbsTol;
opt.RelTol      = param{1}.RelTol;

n_diff          = zeros(n_cells,1);
n_alg           = zeros(n_cells,1);
start_x_index   = 1;
start_xp_index  = 1;

% For each cell, allcoate memory to store external functions used to
% estimate the SOC.
SOC_estimate    = cell(n_cells,1);

% Perform several checks over the cells
for i=1:n_cells
    % Check the daeFormulation flag in case startSimulation is called
    if(param{i}.daeFormulation~=1)
        error(['Make sure that the daeFormulation flag is set to 1 for each cell parameters structure. Cell ', num2str(i),' does not respect this constraint.'])
    end

    % Check if the SOC estimation function handle have been set. In case that
    % the funciton handle has not been defined or it does not have the right
    % number of input arguments, then return empty values.
    if(isempty(param{i}.SOC_estimation_function) || nargin(param{i}.SOC_estimation_function)~=6)
        SOC_estimate{i} = @(a,b,c,d,e,f,g,h,i,j,k)[]; % What exactly is this ?
    else
        SOC_estimate{i} = @SOCestimation;
    end
    param{i}.Nsum      = param{i}.Np + param{i}.Ns + param{i}.Nn;
    param{i}.Nsum_nos  = param{i}.Np + param{i}.Nn;

    % Define the discretization steps.
    param{i}.deltax_al     = 1 / param{i}.Nal;
    param{i}.deltax_p      = 1 / param{i}.Np;
    param{i}.deltax_s      = 1 / param{i}.Ns;
    param{i}.deltax_n      = 1 / param{i}.Nn;
    param{i}.deltax_cu     = 1 / param{i}.Ncu;

    % Compute the indices used to store the positions of the differential
    % and algebraic variables.
    param{i} = computeVariablesIndices(param{i});

    % Preallocate the differentiation matrices used for the solid phase
    % potential. This can be done here because such matrices are considered
    % to be time invariant.
    param{i} = solidPhaseDifferentiationMatrices(param{i});

    % Initialize Param.I_density or Param.P_density, using the value of the current density/power density (set in lines of code a few lines above
    if param{1}.OperatingMode==1 || param{1}.OperatingMode==4
        param{i}.I_density = param{1}.getCurrentDensity(0,t0,tf,param{1}.extraData);
    elseif param{1}.OperatingMode==2 || param{1}.OperatingMode==5
        param{i}.P_density = param{1}.getPowerDensity(0,t0,tf,param{1}.extraData);
    else
        param{i}.I_density=0; % this ia dummy value
        param{i}.P_density =0;
    end

    %% Initial conditions
    % Initial conditions for the differential states
    [cs_average_init, ce_init, T_init, film_init, Q_init, n_diff(i)] =   differentialInitialConditions(param{i});

    % We are now using the 'initialise_model' function for algebraicInitialConditions (refer to code a few lines below)
    n_alg(i)=param{1}.Np+param{1}.Nn+param{1}.Np+param{1}.Nn+param{1}.Nsum+param{1}.Nn+1;

    % Store the number of differential and algebraic variables for each cell.
    param{i}.ndiff = n_diff(i);
    param{i}.nalg  = n_alg(i);

    if((Y0_existence==0) && (YP0_existence==0))
        % Solve the algebraic equations to find a set of semi-consistent initial
        % conditions for the algebraic equations. This will help the DAE solver as a warm startup.
        [init_point, n_alg(i), ~] = initialise_model(param{i}); %init_point is entirely an algebraic vector

        % Build the initial values array for the integrator
        Yt0 = [ce_init;cs_average_init;T_init;film_init;Q_init;init_point];
        Y0  = [Y0;Yt0];
        YP0 = [YP0;zeros(size(Yt0))];
    end

    % The x_index variable will be used in the battery model file
    % for indexing purposes.
    param{i}.x_index    = (start_x_index:n_diff(i)+n_alg(i)+start_x_index-1);
    param{i}.xp_index   = (start_xp_index:n_diff(i)+start_xp_index-1);

    % Update the starting x_index value for the (possible) next cell
    start_x_index       = n_diff(i)+n_alg(i)+start_x_index;
    start_xp_index      = n_diff(i)+n_alg(i)+start_xp_index;
end

if(n_cells==1)
    nc = ' cell';
else
    nc = ' cells';
end

% Empty the used arrays
ce_t            = cell(n_cells,1);
cs_bar_t        = cell(n_cells,1);
T_t             = cell(n_cells,1);
jflux_t         = cell(n_cells,1);
Phis_t          = cell(n_cells,1);
Phie_t          = cell(n_cells,1);
cs_star_t       = cell(n_cells,1);
t_tot           = cell(n_cells,1);
Qrev_t          = cell(n_cells,1);
Qrxn_t          = cell(n_cells,1);
Qohm_t          = cell(n_cells,1);
SOC_t           = cell(n_cells,1);
Voltage_t       = cell(n_cells,1);
SOC_estimated_t = cell(n_cells,1);
film_t          = cell(n_cells,1);
js_t            = cell(n_cells,1);
R_int_t         = cell(n_cells,1);
curr_density_t  = cell(n_cells,1);
Up_t            = cell(n_cells,1);
Un_t            = cell(n_cells,1);
etap_t          = cell(n_cells,1);
etan_t          = cell(n_cells,1);
dudtp_t         = cell(n_cells,1);
dudtn_t         = cell(n_cells,1);
Q_t             = cell(n_cells,1);
yp_original     = YP0';

% This flag is used to notify the reason of the simulation stop. If 0
% everything went well.
exit_reason     = 0;

% Define the structure to be passed to the residual function
ida_user_data.param  = param;
ida_user_data.t0     = t0;
ida_user_data.tf     = tf;

% Define algebraic and differential variables.
% id:1-> differential variables,
% id:0-> algebraic variables.
id = [];
constraint_vector=[];
for i=1:n_cells
    id = [id;ones(n_diff(i),1);zeros(n_alg(i),1)];
    temp_constraint_vector = zeros(n_diff(i)+n_alg(i),1);
    temp_constraint_vector(param{i}.Phis_indices)=1.0; % enforce positivity on solid phase potential in all nodes
    constraint_vector=[constraint_vector;temp_constraint_vector];
end

JacFun = [];

if (param{1}.OperatingMode==4||param{1}.OperatingMode==5)
    param{1}.JacobianFunction=[];
end

% This statement checks if the user wants to make use of the Jacobian
% matrix and (if yes) it has been already provided or not as part of the parameters structure.
if(param{1}.UseJacobian==1 && isempty(param{1}.JacobianFunction))
    % If the user wants to make use of the Jacobian, but it was not
    % provided in the parameters structure, then evaluate a new Jacobian.
    %     disp('Evaluating the analytical form of the Jacobian matrix. Please wait...')

    % Import casadi framework
    import casadi.*
    % Define the symbolic variables.
    xsym    = SX.sym('x',[sum(n_diff)+sum(n_alg),1]);
    xpsym   = SX.sym('xp',[sum(n_diff)+sum(n_alg),1]);
    cj      = SX.sym('cj',1);

    % Get the model equations written in an implicit form in a symbolic way.
    [dx_tot, ~, ~] = batteryModel(0,xsym,xpsym,ida_user_data);

    % Evaluate the Jacobian matrix. (Please refer to the Sundials guide for
    % further information about the Jacobian structure).
    J = jacobian(dx_tot,xsym) + cj*jacobian(dx_tot,xpsym);

    % Define a function for the Jacobian evaluation for a given set of
    % differential and algebraic variables.
    JacFun = Function('fJ',{xsym,cj},{J});

    % Store the function into a structure such that IDA will use it for the
    % evaluation of the Jacobian matrix (see the definition of the function
    % djacfn at the end of this file).
    ida_user_data.fJ = JacFun;

    % Define the options for Sundials
    options = IDASetOptions('RelTol', opt.RelTol,...
        'AbsTol'        , opt.AbsTol,...
        'MaxNumSteps'   , 1500,...
        'VariableTypes' , id,...
        'UserData'      , ida_user_data,...
        'JacobianFn'    , @djacfn,...
        'LinearSolver'  , 'Dense',...
        'ConstraintTypes', constraint_vector);

elseif(param{1}.UseJacobian==1 && ~isempty(param{1}.JacobianFunction))
    % If the Jacobian needs to be used and it has also been provided in the
    % parameters structure, use it directly.

    %     disp('Analytical function of the Jacobian matrix provided by the user.')
    % If the Jacobian has been provided from the user, use it directly.
    JacFun = param{1}.JacobianFunction;

    % Pass this function pointer to the routine that IDA will call for the
    % evaluation of the Jacobian values.
    ida_user_data.fJ = JacFun;

    % Define the options for Sundials
    options = IDASetOptions('RelTol', opt.RelTol,...
        'AbsTol'        , opt.AbsTol,...
        'MaxNumSteps'   , 1500,...
        'VariableTypes' , id,...
        'UserData'      , ida_user_data,...
        'JacobianFn'    , @djacfn,...
        'LinearSolver'  , 'Dense',...
        'ConstraintTypes', constraint_vector);
else
    % In this case the user does not want to make use of the Jacobian
    % matrix. A numerical approximation will be calculated instead.

    % Define the options for Sundials
    options = IDASetOptions('RelTol', opt.RelTol,...
        'AbsTol'        , opt.AbsTol,...
        'MaxNumSteps'   , 1500,...
        'VariableTypes' , id,...
        'UserData'      , ida_user_data,...
        'LinearSolver'  , 'Dense',...
        'ConstraintTypes', constraint_vector);
end

% Initialise solver
IDAInit(@batteryModel,t0,Y0,YP0,options);

if param{i}.Scope == 1
    disp(['Finding a set of consistent ICs for ',num2str(n_cells),nc,' battery pack. Please wait..'])
end

% Find consistent initial conditions
[status, yy, ~] = IDACalcIC(t0+10,'FindAlgebraic'); % might have to change the 10 to a different horizon
% Init the starting integration time
t = t0;

% Store in the results the initial states values.
y = yy';

[ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,Up_t,Un_t,R_int_t,curr_density_t,Voltage_t,SOC_estimated_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,~,dudtp_t, dudtn_t,t_tot] =...
    storeSimulationResults(n_cells,ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,curr_density_t,Voltage_t,SOC_estimated_t,Up_t,Un_t,R_int_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,dudtp_t, dudtn_t, t_tot, y, t,SOC_estimate,t0,tf, param);

sim_time = 0;  % Simulation time (i.e. wall) for reporting purposes only. Not used in controlling solver or time-loop
% Loop until the integration time reaches tf.
while(t<tf)
    %% Check stop conditions for each cell
    for i=1:n_cells
        Phis_pos_cc_t = 1.5*Phis_t{i}(end,1) - 0.5*Phis_t{i}(end,2);
        Phis_neg_cc_t = 1.5*Phis_t{i}(end,end) - 0.5*Phis_t{i}(end,end-1);
        voltage = Phis_pos_cc_t - Phis_neg_cc_t;
        Sout    = internalSOCestimate(cs_bar_t,param,i);
        max_layer_temperature = max(T_t{1}(end,:));
        % Break conditions.
        if(voltage<param{i}.CutoffVoltage)
            if param{i}.suppress_status_prints == 0
                fprintf('\nCell #%d  below its Cutoff voltage. Stopping ...\n',i);
            end
            if exit_reason==0
                exit_reason=[];
            end
            exit_reason = [exit_reason;1];
        end

        if(voltage>param{i}.CutoverVoltage)
            if param{i}.suppress_status_prints == 0
                fprintf('\nCell #%d  above its Cutover voltage. Stopping ...\n',i);
            end
            if exit_reason==0
                exit_reason=[];
            end
            exit_reason = [exit_reason;2];
        end

        if(Sout<param{i}.CutoffSOC)
            if param{i}.suppress_status_prints == 0
                fprintf('\nCell #%d  below its Cutoff SOC. Stopping ...\n',i);
            end
            if exit_reason==0
                exit_reason=[];
            end
            exit_reason = [exit_reason;3];
        end

        if(Sout>param{i}.CutoverSOC)
            if param{i}.suppress_status_prints == 0
                fprintf('\nCell #%d  above its Cutover SOC. Stopping ...\n',i);
            end
            if exit_reason==0
                exit_reason=[];
            end
            exit_reason = [exit_reason;4];
        end

        if(max_layer_temperature > 0.99*param{i}.Tmax)
            if param{i}.suppress_status_prints == 0
                fprintf('\nCell #%d  above its Maximum Permitted Temperature. Stopping ...\n',i);
            end
            if exit_reason==0
                exit_reason=[];
            end
            exit_reason = [exit_reason;5];
        end

        if(param{i}.enable_csneg_Saturation_limit == 1)
            if(param{i}.SolidPhaseDiffusion~=3)
                max_cs_surface_neg = max(cs_bar_t{i}(end,param{i}.Np+1:end));
            else
                max_cs_surface_neg = max(cs_bar_t{i}(end, (param{i}.Np*param{i}.Nr_p) +1:end));
            end

            if (max_cs_surface_neg > param{i}.cs_sat_thresh*param{i}.cs_neg_saturation)
                if param{i}.suppress_status_prints == 0
                    fprintf('\n\nCell #%d  above its safe surface concentration (saturation threshold) Stopping ...\n\n',i);
                end
                if exit_reason==0
                    exit_reason=[];
                end
                exit_reason = [exit_reason;6];
            end
        end
    end

    if(exit_reason~=0)
        break;
    end

    %% Solve the set of DAEs
    % The solver IDA is used to solve the resulting set of DAEs. Please
    % refer to IDA manual for more information about syntax and its usage.
    tic
    [~, t, y]   = IDASolve(tf,'OneStep');
    sim_time    = sim_time+toc;
    y           = y';
    % Store derivative info at each time step
    yp_original = [yp_original;IDAGet('DerivSolution',t,1)'];

    [ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,Up_t,Un_t,R_int_t,curr_density_t,Voltage_t,SOC_estimated_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,tot_voltage,dudtp_t, dudtn_t,t_tot] =...
        storeSimulationResults(n_cells,ce_t,cs_bar_t,T_t,jflux_t,Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t,curr_density_t,Voltage_t,SOC_estimated_t,Up_t,Un_t,R_int_t,Qrev_t,Qrxn_t,Qohm_t,Q_t,dudtp_t, dudtn_t, t_tot, y, t,SOC_estimate,t0,tf, param);

    % If the output scope is active, show additional information to the user
    if(param{1}.Scope==1)
        if(n_cells==1)
            temperature     = T_t{1}(end,end); % If applying a discretised heat model within each layer, this will need to be changed to the max value of temperature within
            % If Fick's law of diffusion is used, before to evaluate the
            % SOC, it is necessary to compute the average solid
            % concentration in each particle.
            Sout = internalSOCestimate(cs_bar_t,param,1);
            clc
            fprintf(['No. of cells in the pack \t',num2str(n_cells),'\n']);
            fprintf(['Time \t\t\t\t\t',num2str(t),' s\n']);
            % If potentiostatic mode is running, applied current comes as
            % solution of DAEs. Otherwise it is provided by the user.
            if(param{1}.OperatingMode==1) || (param{1}.OperatingMode==4)
                fprintf(['Applied current density \t\t',num2str(y(end)),' A/m^2\n']);
            elseif (param{1}.OperatingMode==2) || (param{1}.OperatingMode==5)
                fprintf(['Applied power density \t\t',num2str(param{1}.getPowerDensity(t,t0,tf,param{1}.extraData)),' W/m^2\n']);
            end
            fprintf(['Voltage \t\t\t\t',          num2str((1.5*Phis_t{1}(end,1)-0.5*Phis_t{1}(end,2)) - (1.5*Phis_t{1}(end,end) - 0.5*Phis_t{1}(end,end-1))),   ' V\n']);
            fprintf(['Temperature \t\t\t',        num2str(temperature),                           ' K\n']);
            fprintf(['SOC \t\t\t\t\t',            num2str(Sout),                                  ' %% \n']);
            fprintf(['Cutoff Voltage \t\t\t',     num2str(param{1}.CutoffVoltage),                ' V\n']);
            fprintf(['Cutover Voltage \t\t',      num2str(param{1}.CutoverVoltage),               ' V\n']);
            fprintf(['Internal Resistance \t',    num2str(R_int_t{1}(end)),                       ' Ohm m^2\n']);
            fprintf(['Absolute tolerance \t\t',   num2str(param{1}.AbsTol),                       '\n']);
            fprintf(['Relative tolerance \t\t',   num2str(param{1}.RelTol),                       '\n']);
            fprintf(['Initial int. time \t\t',    num2str(t0),                                    ' s\n']);
            fprintf(['Final int. time \t\t',      num2str(tf),                                    ' s\n']);
            fprintf(['N. of unknowns \t\t\t',     num2str(length(y)),                             ' \n']);
        else
            clc
            fprintf(['No. of cells in the pack \t',num2str(n_cells),'\n']);
            fprintf(['Time \t\t\t\t\t',num2str(t),' s\n']);
            % If potentiostatic mode is running, applied current comes as
            % solution of DAEs. Otherwise it is provided by the user.
            if(param{1}.OperatingMode==1) || (param{1}.OperatingMode==4)
                fprintf(['Applied current density \t\t',num2str(y(end)),' A/m^2\n']);
            elseif (param{1}.OperatingMode==2) || (param{1}.OperatingMode==5)
                fprintf(['Applied power density \t\t',num2str(param{1}.getPowerDensity(t,t0,tf,param{1}.extraData)),' W/m^2\n']);
            end

            fprintf(['Voltage \t\t\t\t',          num2str(tot_voltage),       ' V\n']);
            fprintf(['Absolute tolerance \t\t',   num2str(param{1}.AbsTol),   ' \n']);
            fprintf(['Relative tolerance \t\t',   num2str(param{1}.RelTol),   ' \n']);
            fprintf(['Initial int. time \t\t',    num2str(t0),                ' s\n']);
            fprintf(['Final int. time \t\t',      num2str(tf),                ' s\n']);
            fprintf(['N. of unknowns \t\t\t',     num2str(length(y)),         ' \n']);
        end
    end
end

% disp(['Elasped time: ',num2str(sim_time),' s']);

% Interpolate for fixed time step values
t_tot_original = t_tot;

% Build the time vector used for interpolation
time_vector = (t0:param{i}.sim_datalog_interval:tf);

% In case of the simulation has stopped before the final time set by the
% user, change the tf variable in order to interpolate only available values.

if(t<tf)
    % Set final time equal to the last integration step
    tf          = t;
    % Redefine the time vector used for interpolation
    time_vector = (t0:param{i}.sim_datalog_interval:tf);
end

% If at least one integration step has been done, retreive the first order
% time derivative information. Otherwise use the initial data.
if(time_vector(end)>t0)
    % Retreive derivative information at the last time step
    yp          = interp1(t_tot{i},yp_original,time_vector(end))';
    % After interpolating, delete all the other data.
    yp_original = yp_original(end,:)';
else
    % If the integration step carried out by SUNDIALS is less than the
    % parametrized step size, then return the initial data as set of initial states.
    yp          = YP0;
    yp_original = YP0;
end

% Free memory allocated by IDA solver
IDAFree

% These variables will be used to store the original results of the integration process.
Phis_t_o            = cell(n_cells,1);
Phie_t_o            = cell(n_cells,1);
ce_t_o              = cell(n_cells,1);
cs_star_t_o         = cell(n_cells,1);
cs_average_t_o      = cell(n_cells,1);
jflux_t_o           = cell(n_cells,1);
SOC_t_o             = cell(n_cells,1);
T_t_o               = cell(n_cells,1);
Voltage_t_o         = cell(n_cells,1);
SOC_estimated_t_o   = cell(n_cells,1);
film_t_o            = cell(n_cells,1);
js_t_o              = cell(n_cells,1);
R_int_t_o           = cell(n_cells,1);
Up_t_o              = cell(n_cells,1);
Un_t_o              = cell(n_cells,1);
dudtp_t_o           = cell(n_cells,1);
dudtn_t_o           = cell(n_cells,1);
% Qrev_t_o            = cell(n_cells,1);
% Qrxn_t_o            = cell(n_cells,1);
% Qohm_t_o            = cell(n_cells,1);
etap_t_o            = cell(n_cells,1);
etan_t_o            = cell(n_cells,1);
app_current_t_o     = cell(n_cells,1);
Q_t_o               = cell(n_cells,1);
y                   = [];
y_original          = [];
for i=1:n_cells
    % Save the overpotentials
    etap_t{i} = Phis_t{i}(:,1:param{i}.Np)-Phie_t{i}(:,1:param{i}.Np)-Up_t{i};
    if(param{i}.EnableAgeing==1)
        etan_t{i} = Phis_t{i}(:,param{i}.Np+1:end)-Phie_t{i}(:,param{i}.Np+param{i}.Ns+1:end)-Un_t{i} - param{1}.F*jflux_t{i}(:,param{i}.Np+1:end).*(param{i}.R_SEI+film_t{i}./param{i}.k_n_aging);
    else
        etan_t{i} = Phis_t{i}(:,param{i}.Np+1:end)-Phie_t{i}(:,param{i}.Np+param{i}.Ns+1:end)-Un_t{i};
    end
    if(param{i}.sim_datalog_interval>0)
        % Store original results
        Phis_t_o{i}            = Phis_t{i};
        Phie_t_o{i}            = Phie_t{i};
        ce_t_o{i}              = ce_t{i};
        cs_star_t_o{i}         = cs_star_t{i};
        cs_average_t_o{i}      = cs_bar_t{i};
        jflux_t_o{i}           = jflux_t{i};
        SOC_t_o{i}             = SOC_t{i};
        T_t_o{i}               = T_t{i};
        Voltage_t_o{i}         = Voltage_t{i};
        SOC_estimated_t_o{i}   = SOC_estimated_t{i};
        film_t_o{i}            = film_t{i};
        js_t_o{i}              = js_t{i};
        R_int_t_o{i}           = R_int_t{i};
        app_current_t_o{i}     = curr_density_t{i};
        Up_t_o{i}              = Up_t{i};
        Un_t_o{i}              = Un_t{i};
        dudtp_t_o{i}           = dudtp_t{i};
        dudtn_t_o{i}           = dudtn_t{i};
        etap_t_o{i}            = etap_t{i};
        etan_t_o{i}            = etan_t{i};
        % Qrev_t_o{i}            = Qrev_t{i};
        % Qrxn_t_o{i}            = Qrxn_t{i};
        % Qohm_t_o{i}            = Qohm_t{i};
        Q_t_o{i}               = Q_t{i};

        if(time_vector(end)>t0)
            % Interpolate the results
            Phis_t{i}          = interp1(t_tot{i},Phis_t{i},time_vector');
            Phie_t{i}          = interp1(t_tot{i},Phie_t{i},time_vector');
            ce_t{i}            = interp1(t_tot{i},ce_t{i},time_vector');
            cs_star_t{i}       = interp1(t_tot{i},cs_star_t{i},time_vector');
            cs_bar_t{i}        = interp1(t_tot{i},cs_bar_t{i},time_vector');
            jflux_t{i}         = interp1(t_tot{i},jflux_t{i},time_vector');
            SOC_t{i}           = interp1(t_tot{i},SOC_t{i},time_vector');
            SOC_estimated_t{i} = interp1(t_tot{i},SOC_estimated_t{i},time_vector');
            Voltage_t{i}       = interp1(t_tot{i},Voltage_t{i},time_vector');
            film_t{i}          = interp1(t_tot{i},film_t{i},time_vector');
            js_t{i}            = interp1(t_tot{i},js_t{i},time_vector');
            R_int_t{i}         = interp1(t_tot{i},R_int_t{i},time_vector');
            T_t{i}             = interp1(t_tot{i},T_t{i},time_vector');
            curr_density_t{i}  = interp1(t_tot{i},curr_density_t{i},time_vector');
            Up_t{i}            = interp1(t_tot{i},Up_t{i},time_vector');
            Un_t{i}            = interp1(t_tot{i},Un_t{i},time_vector');
            % Qrev_t{i}          = interp1(t_tot{i},Qrev_t{i},time_vector');
            % Qrxn_t{i}          = interp1(t_tot{i},Qrxn_t{i},time_vector');
            % Qohm_t{i}          = interp1(t_tot{i},Qohm_t{i},time_vector');
            etap_t{i}          = interp1(t_tot{i},etap_t{i},time_vector');
            etan_t{i}          = interp1(t_tot{i},etan_t{i},time_vector');
            dudtp_t{i}         = interp1(t_tot{i},dudtp_t{i},time_vector');
            dudtn_t{i}         = interp1(t_tot{i},dudtn_t{i},time_vector');
            Q_t{i}             = interp1(t_tot{i},Q_t{i},time_vector');
            t_tot{i}           = time_vector';
        end
    end
    % Store results. If integration steps are enabled, store the interpolated data.
    results.Phis{i}                        = Phis_t{i};
    results.Phie{i}                        = Phie_t{i};
    results.ce{i}                          = ce_t{i};
    results.cs_surface{i}                  = cs_star_t{i};
    results.cs_average{i}                  = cs_bar_t{i};
    results.time{i}                        = t_tot{i};
    results.int_internal_time{i}           = t_tot_original{i};
    results.ionic_flux{i}                  = jflux_t{i};
    results.side_reaction_flux{i}          = js_t{i};
    results.SOC{i}                         = SOC_t{i};
    results.SOC_estimated{i}               = SOC_estimated_t{i};
    results.Voltage{i}                     = Voltage_t{i};
    results.Temperature{i}                 = T_t{i};
    results.Qrev{i}                        = Qrev_t{i};
    results.Qrxn{i}                        = Qrxn_t{i};
    results.Qohm{i}                        = Qohm_t{i};
    results.film{i}                        = film_t{i};
    results.R_int{i}                       = R_int_t{i};
    results.Up{i}                          = Up_t{i};
    results.Un{i}                          = Un_t{i};
    results.etap{i}                        = etap_t{i};
    results.etan{i}                        = etan_t{i};
    results.dudtp{i}                       = dudtp_t{i};
    results.dudtn{i}                       = dudtn_t{i};
    results.Q{i}                           = Q_t{i};
    results.parameters{i}                  = param{i};
    results.JacobianFun                    = JacFun;

    % Store original data.
    results.original.Phis{i}               = Phis_t_o{i};
    results.original.Phie{i}               = Phie_t_o{i};
    results.original.ce{i}                 = ce_t_o{i};
    results.original.cs_surface{i}         = cs_star_t_o{i};
    results.original.cs_average{i}         = cs_average_t_o{i};
    results.original.ionic_flux{i}         = jflux_t_o{i};
    results.original.side_reaction_flux{i} = js_t_o{i};
    results.original.SOC{i}                = SOC_t_o{i};
    results.original.SOC_estimated{i}      = SOC_estimated_t_o{i};
    results.original.Voltage{i}            = Voltage_t_o{i};
    results.original.Temperature{i}        = T_t_o{i};
    results.original.film{i}               = film_t_o{i};
    results.original.R_int{i}              = R_int_t_o{i};
    results.original.Up{i}                 = Up_t_o{i};
    results.original.Un{i}                 = Un_t_o{i};
    results.original.etap{i}               = etap_t_o{i};
    results.original.etan{i}               = etan_t_o{i};
    results.original.Q{i}                  = Q_t_o{i};
    results.original.parameters{i}         = param_original{i};

    % Store initial states data
    y           = [y;ce_t{i}(end,:)';cs_bar_t{i}(end,:)';T_t{i}(end,:)';film_t{i}(end,:)';Q_t{i}(end,:)';jflux_t{i}(end,:)';Phis_t{i}(end,:)';Phie_t{i}(end,:)';js_t{i}(end,:)';curr_density_t{i}(end)];
    % Store initial states original data
    y_original  = [y_original;ce_t_o{i}(end,:)';cs_average_t_o{i}(end,:)';T_t_o{i}(end,:)';film_t_o{i}(end,:)';Q_t_o{i}(end,:)';jflux_t_o{i}(end,:)';Phis_t_o{i}(end,:)';Phie_t_o{i}(end,:)';js_t_o{i}(end,:)';app_current_t_o{i}(end)];
end

% Store the array of last results
results.Y                           = y;
results.YP                          = yp;

results.original.Y                  = y_original;
results.original.YP                 = yp_original;

results.original.initialState.Y     = y_original;
results.original.initialState.YP    = yp_original;

results.initialState.Y              = y;
results.initialState.YP             = yp;

% Store simulation time
results.simulation_time             = sim_time;

% Exit reason
results.exit_reason                 = exit_reason;

% Check the operating mode and store the results accordingly.
if(param{1}.OperatingMode==2)
    results.power_density          = param{1}.P_density * ones(size(t_tot{1},1),1);
    results.original.power_density = param{1}.P_density * ones(size(t_tot_original{1},1),1);
    % Variable current profile
elseif(param{1}.OperatingMode==5)
    results.power_density          = param{1}.getPowerDensity(t_tot{1},t0,tf,param{1}.extraData);
    results.original.power_density = param{1}.getPowerDensity(t_tot_original{1},t0,tf,param{1}.extraData);
end

results.curr_density          = curr_density_t{1};
results.original.curr_density = app_current_t_o{1};
end

function estimate = SOCestimation(t,t0,tf,param,ce_t,cs_bar_t,cs_star_t,Phie_t,Phis_t,jflux_t,T_t)
% Build the states struct which will be passed to the function
states.ce           = ce_t(end,:);
states.cs_average   = cs_bar_t(end,:);
states.cs_surface   = cs_star_t(end,:);
states.Phie         = Phie_t(end,:);
states.Phis         = Phis_t(end,:);
states.ionic_flux   = jflux_t(end,:);
states.Temperature  = T_t(end,:);
% Call the estimation procedure
estimate = param.SOC_estimation_function(t,t0,tf,states,param.extraData,param);
end

% This function is used to get a measurement of the SOC according to the
% internal states. This function assumes that all the states are measurable.
function Sout = internalSOCestimate(cs_average_t,param,i)

% Check if Fick's law of diffusion is used. This is required to define the correct way how to evaluate the SOC.
% by using cs_average in neg electrode
if(param{i}.SolidPhaseDiffusion~=3)
    cs_average = sum(cs_average_t{i}(end,param{i}.Np+1:end))/(param{i}.Nn);
else
    cs_average = sum(cs_average_t{i}(end, (param{i}.Np*param{i}.Nr_p) +1:end))/(param{i}.Nn*param{i}.Nr_n);
end
Sout = 100*((cs_average/param{i}.cs_maxn) - param{i}.theta_min_neg)/(param{i}.theta_max_neg- param{i}.theta_min_neg); % cell's soc in percentage
end

% This function is used to evaluate the Jacobian Matrix of the P2D model.
function [J, flag, new_data] = djacfn(t, y, yp, rr, cj, data)

% Extract the function object (representing the Jacobian)
fJ    = data.fJ;

% Evaluate the Jacobian with respect to the present values of the states and their time derivatives.
try
    J = full(fJ(y,cj));
catch
    J = [];
end

% Return dummy values
flag        = 0;
new_data    = [];
end