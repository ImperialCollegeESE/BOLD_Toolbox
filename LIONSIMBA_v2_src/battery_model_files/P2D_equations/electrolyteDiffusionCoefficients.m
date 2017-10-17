function Deff = electrolyteDiffusionCoefficients(ce,T,param,batterySection)
% ELECTROLYTEDIFFUSIONCOEFFICIENTS  Evaluates the diffusion coefficients for the electrolyte phase [m^2/s].
%
%   [Deff_p, Deff_s, Deff_n] = ELECTROLYTEDIFFUSIONCOEFFICIENTS(ce,T,param) evaluates
%   the diffusion coefficients for the anode, separator and cathode of the
%   battery.
%
%   Note that this is an interface for the main program. The authors
%   suggest to maintain the name of the function and its signature, while
%   modifying only the body of the script.
%
%   The diffusion coefficients can be evaluated in isothermal case
%   (param.TemperatureEnabled=0) or adiabatic case
%   (param.TemperatureEnabled=1).
%
%   You can modify the way that the diffusion coefficients are computed, as
%   function of electrolyte concentration and temperature. The main script
%   will pass also the param array.
%
%   The user can modify this script to meet particular requirements.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.]
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(param.TemperatureEnabled==1)
    switch(batterySection)
        case 'p'
            Deff = param.eps_p^param.brugg_p*1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
        case 's'
            Deff = param.eps_s^param.brugg_s*1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
        case 'n'
            Deff = param.eps_n^param.brugg_n*1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
    end
else
    switch(batterySection)
        case 'p'
            Deff = repmat(param.Dp*param.eps_p^param.brugg_p,param.Np,1);
        case 's'
            Deff = repmat(param.Ds*param.eps_s^param.brugg_s,param.Ns,1);
        case 'n'
            Deff = repmat(param.Dn*param.eps_n^param.brugg_n,param.Nn,1);
    end
end

end
