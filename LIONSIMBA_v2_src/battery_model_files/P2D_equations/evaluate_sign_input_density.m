function sign_input_density = evaluate_sign_input_density(param)
% returns the sign of input current/power density as per operating mode

if param.OperatingMode==1 || param.OperatingMode==4
    sign_input_density = sign(param.I_density);
elseif param.OperatingMode==2 || param.OperatingMode==5
    sign_input_density = sign(param.P_density);
elseif param.OperatingMode==3
    sign_input_density = 1; % dummy value for CV mode (since usually we have CV charging)
else
    error('Not a valid operating mode.');
end

end
