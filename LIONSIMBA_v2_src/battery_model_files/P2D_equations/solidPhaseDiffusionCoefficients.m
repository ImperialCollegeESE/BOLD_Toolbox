function [Dps_eff, Dns_eff] = solidPhaseDiffusionCoefficients(T,param)
% SOLIDPHASEDIFFUSIONCOEFFICIENTS evaluates diffusion coefficients of the solid phase [m^2 /s].
% The user may modify the script to meet specific requirements.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(param.TemperatureEnabled==1)
    Dps_eff     = param.Dps*exp(-param.EaDps/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
else
    Dps_eff     = param.Dps*ones(param.Np,1);
end

if(param.TemperatureEnabled==1)
    Dns_eff     = param.Dns*exp(-param.EaDns/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
else
    Dns_eff     = param.Dns*ones(param.Nn,1);
end

end
