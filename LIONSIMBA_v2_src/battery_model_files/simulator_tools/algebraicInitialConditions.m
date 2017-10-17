function [x0_alg, n_alg] =   algebraicInitialConditions(param)
% algebraicInitialConditions returns the number of algebraic variables
% DEPRECATED in favor of initialise_model.m

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess for the algebraic variables
jflux_init          = [-0.43e-5*ones(param.Np,1);0.483e-5*ones(param.Nn,1)];
Phis_init           = [4.2*ones(param.Np,1);0.074*ones(param.Nn,1)];
Phie_init           = zeros(param.Nsum,1);
js_init             = 0.483e-5*ones(param.Nn,1);

% Never mind: these variables returned are not used, since we have analytical initialisation!
% DEPRECATED in favor of initialise_model.m (analytical initialisation)
if param.OperatingMode==1 || param.OperatingMode==4
    I_density       = param.I_density;
elseif param.OperatingMode==2 || param.OperatingMode==5
    I_density       = param.P_density/(Phis_init(1)-Phis_init(end)); % no need for linear interpolation since starting from equilibrium
else
    I_density=0; % dummy value for CV mode
end

% Build the array of algebraic initial conditions
x0_alg              = [
    jflux_init;...
    Phis_init;...
    Phie_init;...
    js_init;...
    I_density
    ];

n_alg = length(x0_alg);
end
