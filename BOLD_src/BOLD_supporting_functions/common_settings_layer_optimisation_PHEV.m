% Authors        : Krishnakumar Gopalakrishnan, Ian D. Campbell, Imperial College London
%                : Davide M. Raimondo, University of Pavia
% copyright year : 2017
% Last Updated   : Tue Oct 17 21 : 44 : 42 CEST 2017
% Licensed       : MIT License
OperatingMode      = 2;    % LIONSIMBA operating mode (constant power density)
cutoffSOC          = 30;   % Cutoff SOC  [%]   If changing SOC limits, remember to change the initial & finishing SOCs for fast charging, if required!
CutoverSOC         = 90;   % Cutover SOC [%]
CutoffVoltage      = 3.50; % Not operate below this value [V]
CutoverVoltage     = 4.20; % Not operate below above value [V]
solidDiffusionMode = 3;    % Enable Fick's law
TemperatureEnabled = 1;    % Flag to turn on/off thermal mode

axial_nodes  = 40; % Number of axial nodes i.e. Np, Nn, Ns along the through-thickness direction
radial_nodes = 15; % Number of radial nodes i.e. Nr_p, Nr_n in the spherical particles
% ref: A.M. Bizeray, S. Zhao, S.R. Duncan, D.A. Howey, Lithium-ion battery thermal-electrochemical model-based state estimation using orthogonal collocation and a modified extended Kalman filter, Journal of Power Sources, Volume 296, 20 November 2015, Pages 400-412, ISSN 0378-7753, https://doi.org/10.1016/j.jpowsour.2015.07.019.
% (http://www.sciencedirect.com/science/article/pii/S0378775315300677)

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: