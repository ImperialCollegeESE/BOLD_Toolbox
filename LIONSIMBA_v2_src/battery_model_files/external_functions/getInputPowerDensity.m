function P_density = getInputPowerDensity(t,t0,tf,extra)
% getInputPowerDensity returns the value of the input current density as a function of time.
%
%       P = GETINPUTPOWERDENSITY(t,t0,tf,extra)
%
%       Inputs:
%               - t     : value of the current time step
%               - t0    : initial integration time
%               - tf    : final integration time
%               - extra : extra parameters
%       Outputs:
%               - P_density     : Applied power density [W/m^2]

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%			        Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define your own linear/nonlinear function of time for the applied power density

% P = (t-t0)/(tf-t0) *(-30) + 0;
% P_density = -90*sin(t+100/100);

% P_density = -10*sin(omega*t+phi);
frequency = 0.1;
P_density = -300*sin((2*pi*frequency)*t+(1/100));
if P_density == 0
    P_density = 1e-2;
end
    
end
