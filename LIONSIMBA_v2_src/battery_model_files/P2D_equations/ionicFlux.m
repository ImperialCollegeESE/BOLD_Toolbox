function [jflux,U_p,U_n,dudt_p,dudt_n,J_s] = ionicFlux(ce,cs_star,Phis,Phie,T,solverFlux,film,param,sign_input_density,I_density)
% IONICFLUX Computes the molar flux density of Li-ions at the electrode-electrolyte interface [mol/(m^2*s)].

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Positive electrode
% Compute the OCV for the positive and negative electrodes.
% [U_p,dudt_p,U_n,dudt_n] = param.OpenCircuitPotentialFunction(cs_star,T,param);
[U_p,dudt_p,U_n,dudt_n] = openCircuitPotential(cs_star,T,param,sign_input_density);
% Compute the reaction rates.
[k_pT, k_nT] = param.ReactionRatesFunction(T,param);

% Positive electrode ion flux
deltap = ((0.5*param.F)./(param.R*T(param.Nal+1:param.Nal+param.Np))).*(Phis(1:param.Np)-Phie(1:param.Np)-U_p);
ip = 2*k_pT.*sqrt(ce(1:param.Np)).*sqrt(cs_star(1:param.Np)).*sqrt(param.cs_maxp-cs_star(1:param.Np));
jnp_calc = ip.* sinh(deltap);

%% Negative electrode

% If ageing is enabled, take into account the SEI resistance
if(param.EnableAgeing==1)
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n -param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging)));
else
    eta_n   = (Phis(param.Np+1:end)-Phie(param.Np+param.Ns+1:end)-U_n);
end

deltan      = ((0.5*param.F)./(param.R*T(param.Nal+param.Np+param.Ns+1:param.Nal+param.Np+param.Ns+param.Nn))).*eta_n;
in          = 2*k_nT.*sqrt(ce(param.Np+param.Ns+1:end)).*sqrt(cs_star(param.Np+1:end)).*sqrt(param.cs_maxn-cs_star(param.Np+1:end));
jnn_calc    = in.* sinh(deltan);

J_s = zeros(param.Nn,1);

% Switch cases when the applied current density is a symbolic variable
if(isa(I_density,'casadi.SX')||isa(I_density,'casadi.SX') && param.EnableAgeing == 1)
    eta_s = Phis(param.Np+1:end) - Phie(param.Np+param.Ns+1:end) - param.Uref_s - param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging));
    % Tafel equation for the side reaction flux.
    alpha   = 0.5*param.F./(param.R*T(param.Nal+param.Np+param.Ns+1:end-param.Ncu));

    if sign_input_density>=0
        J_s=-param.i_0_jside.*(I_density/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F;
    else
        J_s=zeros(param.Nn,1);
    end

    % Through this if_else CASADI statement, it is possible to represent dynamics that switch according to the value of the symbolic quantity current_density
    % J_s = if_else(sign_input_density>=0,-param.i_0_jside.*(I_density/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F,zeros(param.Nn,1));
elseif(param.EnableAgeing == 1 && sign_input_density > 0) % check with Marcello if > or >=0
    eta_s = Phis(param.Np+1:end) - Phie(param.Np+param.Ns+1:end) - param.Uref_s - param.F*solverFlux(param.Np+1:end).*(param.R_SEI+film./(param.k_n_aging));
    % Tafel equation for the side reaction flux.
    alpha   = 0.5*param.F./(param.R*T(param.Nal+param.Np+param.Ns+1:end-param.Ncu));
    J_s = -param.i_0_jside.*(I_density/param.I1C)^param.w.*(exp(-alpha.*eta_s))./param.F;
end
%% Return value
jflux = [jnp_calc;jnn_calc];
