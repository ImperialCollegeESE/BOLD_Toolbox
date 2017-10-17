function [combined_electrode_thickness_per_layer,len_p,len_n] = compute_domain_thicknesses(num_layers,neg_to_pos_len_ratio,param)
% compute_domain_thicknesses computes the individual thicknesses of each domain for a given len/pos ratio
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

    if iscell(param)
        local_param = param{1};
    else
        local_param = param;
    end

    % Outermost electrodes are either a) both Cu electrodes or b) Cu/Al combo for the formulae used here to be valid
    combined_electrode_thickness_per_layer = (local_param.t_stack ...
        - (ceil(0.5*(num_layers + 1))*local_param.len_cu) ...
        - (ceil(0.5*num_layers)*local_param.len_al))/num_layers - local_param.len_s;  % computes (L_pos + L_neg) for one layer

    len_p = combined_electrode_thickness_per_layer/(neg_to_pos_len_ratio + 1);
    len_n = combined_electrode_thickness_per_layer - len_p;

end

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: