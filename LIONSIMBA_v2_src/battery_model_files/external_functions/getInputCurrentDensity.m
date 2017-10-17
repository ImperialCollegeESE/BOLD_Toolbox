function I_density = getInputCurrentDensity(t,t0,tf,extra)
% getInputCurrentDensity  returns the value of the input current density as a function of time.
%
%       I_density = GETINPUTCURRENTDENSITY(t,t0,tf,extra)
%
%       Inputs:
%              t     : value of the current time step
%              t0    : initial integration time
%              tf    : final integration time
%              extra : extra parameters
%       Outputs:
%              I_density : Applied current density [A/m^2]

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%			        Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define your own linear/nonlinear function of time for the applied current density

% I = (t-t0)/(tf-t0) *(-30) + 0;
I_density = -30*sin(t/100);
end
