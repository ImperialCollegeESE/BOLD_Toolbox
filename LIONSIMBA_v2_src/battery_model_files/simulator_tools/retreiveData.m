function [ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t, cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot] = retreiveData(ce_t, cs_barrato_t, T_t, jflux_t, Phis_t, Phie_t,cs_star_t, SOC_t, film_t, js_t, Q_t,t_tot, y, t, param)
% RETREIVEDATA retreives the data which is returned in the results structure.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract differential variables after the integration process
ce_t            = [ce_t;y(param.ce_indices)];
cs_barrato_t    = [cs_barrato_t;y(param.cs_average_indices)];
T_t             = [T_t;y(param.T_indices)];
film_t          = [film_t;y(param.film_indices)];
Q_t             = [Q_t;y(param.Q_indices)];
% Extract the algebraic variables after the integration process
jflux_t         = [jflux_t;y(param.jflux_indices)];
Phis_t          = [Phis_t;y(param.Phis_indices)];
Phie_t          = [Phie_t;y(param.Phie_indices)];
js_t            = [js_t;y(param.js_indices)];

if(param.SolidPhaseDiffusion~=3)
    cs_average = sum(cs_barrato_t(end,param.Np+1:end))/(param.Nn);  % cs_average in neg electrode
else
    cs_average = sum(cs_barrato_t(end, (param.Np*param.Nr_p) +1:end))/(param.Nn*param.Nr_n); % cs_average in neg electrode
end

Sout = 100*((cs_average/param.cs_maxn) - param.theta_min_neg)/(param.theta_max_neg- param.theta_min_neg); % cell-soc percent

cs_star_t       = [cs_star_t;surfaceConcentration(cs_barrato_t(end,:)',jflux_t(end,:)',Q_t(end,:)',T_t(end,:)',param)'];
SOC_t           = [SOC_t;Sout];
t_tot           = [t_tot;t];
end
