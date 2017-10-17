function n_max_numerical_limit = compute_max_layers_numerical_limit(neg_to_pos_len_ratio,param)
% compute_max_layers_numerical_limit returns the possible max value of number of layers that can be used.
% Maintains the same len/pos ratio for all layers. Maintains the same stack thickness too.
% Also ensures that combined domain thicknesses shall not become negative.

% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

    t_stack = param{1}.t_stack;
    len_cu = param{1}.len_cu;
    len_al = param{1}.len_al;
    len_s  = param{1}.len_s;

    n_max_even = floor(2*(t_stack - len_cu)/(len_al + len_cu + 2*len_s));
    n_max_odd  = floor((2*t_stack - len_al - len_cu)/(len_al + len_cu + 2*len_s));

    n_max_numerical_limit = max(n_max_even,n_max_odd);

end
% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t:
