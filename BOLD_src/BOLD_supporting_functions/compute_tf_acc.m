function [vf, tf] = compute_tf_acc(acc_test_choice, acc_time_manuf, vf_manuf)
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

vf_sae       = 20;  % [mph] (SAE J1666 - final speed after acceleration)
acc_time_sae = 6.0; % [sec] (SAE J1666 - final speed after acceleration)

vf_choices       = [vf_manuf,vf_sae]*0.44704;     % [m/s]
acc_time_choices = [acc_time_manuf,acc_time_sae]; % [m/s]
acc_rate_choices = vf_choices./acc_time_choices;  % [m/s^2]
if strcmp(acc_test_choice,'sae_derived')
    [~, min_idx] = min(acc_rate_choices);     % [m/s^2]
    vf           = vf_choices(min_idx);       % [m/s]
    tf           = acc_time_choices(min_idx); % [sec]
elseif strcmp(acc_test_choice,'worst_case')
    [~, max_idx] = max(acc_rate_choices);     % [m/s^2]
    vf           = vf_choices(max_idx);       % [m/s]
    tf           = acc_time_choices(max_idx); % [sec]
else
    error('\nInvalid test choice....\n');
end

end     % End of function

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: