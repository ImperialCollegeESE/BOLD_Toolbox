function [res_dT, rhsT] = thermalModel(~,cs_star_avg,Phis,~,~,T,dT,Up,Un,dUdT_p,dUdT_n,param,I_density)
% THERMALMODEL evaluates the set of equations for (lumped) thermal dynamics.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uocp_avg_pos = sum(Up,1)/param.Np;
% Uocp_avg_neg = sum(Un,1)/param.Nn;
% Uocp_avg = Uocp_avg_pos - Uocp_avg_neg;

[Up_avg,dUdT_p_avg,Un_avg,dUdT_n_avg] = openCircuitPotential(cs_star_avg,T,param,1); %  the final argument is a dummy argument (currently under debug and may be removed in the future)

Up_pos_cc = 1.5*Up_avg(1) - 0.5*Up_avg(2);
Un_neg_cc = 1.5*Un_avg(end) - 0.5*Un_avg(end-1);
cell_avg_OCP = Up_pos_cc - Un_neg_cc;

Phis_pos_cc = 1.5*Phis(1) - 0.5*Phis(2);
Phis_neg_cc = 1.5*Phis(end) - 0.5*Phis(end-1);
V_cell = Phis_pos_cc - Phis_neg_cc;

cell_current = I_density*param.overall_surface_area_for_given_layers;
Q_polarisation = abs(cell_current)*abs(cell_avg_OCP - V_cell);	% shall always be positive

if param.lumped_thermal_version == 1 % No entropic heat generation
    mCpdTdt = -(param.h_lumped*param.tab_area)*(T - param.Tref)...
              + Q_polarisation;

elseif param.lumped_thermal_version == 2 % With entropic heat generation
    dUdT_p_avg = sum(dUdT_p)/param.Np;
    dUdT_n_avg = sum(dUdT_n)/param.Nn;

    sign_input_density = evaluate_sign_input_density(param); % evaluates sign of input current/power density as per operating mode

    if sign_input_density==1
        Q_rev = cell_current*T*(dUdT_p_avg-dUdT_n_avg); % Charging Entropic heat term
    else
        Q_rev = -cell_current*T*(dUdT_p_avg-dUdT_n_avg);  % Discharging Entropic heat term
    end

    mCpdTdt = -(param.h_lumped*param.tab_area)*(T - param.Tref)...
              + Q_polarisation + Q_rev;
else
    error("Incorrect value for choice of thermal model. Exiting..\n")
end

rhsT = mCpdTdt/(param.mass_cell*param.Cp_avg);
% rhsT = 0;  % REMEMBER TO COMMENT THIS OUT;
res_dT  = (dT - rhsT);

end
