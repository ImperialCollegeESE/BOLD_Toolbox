function [k_pT, k_nT] = reactionRates(T,param)
% REACTIONRATES  Reaction rates (k) of cathode and anode [m^2.5/(m^0.5 s)].
% The user may modify this script to meet specific requirements.

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
    k_pT     = param.k_p*exp(-param.Eakip/param.R*(1./T(param.Nal+1:param.Nal+param.Np)-1/param.Tref));
else
    k_pT     = param.k_p;
end

if(param.TemperatureEnabled==1)
    k_nT     = param.k_n*exp(-param.Eakin/param.R*(1./T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn)-1/param.Tref));
else
    k_nT     = param.k_n;
end

end
