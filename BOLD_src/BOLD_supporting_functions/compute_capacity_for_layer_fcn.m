function neg_electrode_capacity_Ah = compute_capacity_for_layer_fcn(num_layers,param_100pct_soc,param_0pct_soc,len_n)
% returns the capacity of both electrodes for given layer choice
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

ne = 1; % No. of electrons transferred in reaction (for Li-ion, ne = 1)
neg_electrode_volume      = num_layers*param_100pct_soc.surface_area_per_face_Northrop_cell*len_n;
neg_electrode_capacity_Ah = (param_100pct_soc.cs_n_init - param_0pct_soc.cs_n_init)*neg_electrode_volume*param_100pct_soc.F*ne/3600;
neg_electrode_capacity_Ah = neg_electrode_capacity_Ah*(1-param_100pct_soc.eps_n - param_100pct_soc.eps_fi(3));

end

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: