function P_wheels_cruise = compute_cruise_power(vehicle_speed, Mv, Cd, Av)
% Returns a vector of vehicle's power demand (watts) for a given speed (vs time) profile without considering power required to accelerate its mass.
% Positive power represents propulsion/traction. Negative power represents braking.
% Note: this is opposite to LIONSIMBA conventions for discharge/charge
% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License

Cr  = 0.01;   % Rolling resistance coefficient
Z   = 0.0;    % road grade, assumed flat
g   = 9.81;   % m/s^2
rho = 1.2041; % kg/m^3, At 20°C and 101.325 kPa for dry air (https://en.wikipedia.org/wiki/Density_of_air)

P_drag          = 0.5*rho*Cd*Av*vehicle_speed.^3;   % power required to overcome air resistance
P_friction      = Cr*Mv*g*vehicle_speed;            % power required to overcome friction
P_gradient      = Mv*g*Z*vehicle_speed;             % Incline (i.e. road gradient)
P_wheels_cruise = P_drag + P_friction + P_gradient; % Total power required without acceleration power

end

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: