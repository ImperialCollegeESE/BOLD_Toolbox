function checkEnvironment(param,num_arg)
% CHECKENVIRONMENT checks if the required tools are present in the host system.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check for input arguments
if(num_arg==0)
    error('Type help startSimulation')
end

if(num_arg == 3 && param{1}.OperatingMode == 1)
    error('For constant input currents, please provide the current density as a parameter to startSimulation')
end

% Check if SUNDIALS is installed correctly
% disp('Checking for SUNDIALS availability...');
result = exist('IDAFree');
if(result~=2)
    error('The SUNDIALS software does not seem to be present in the search path or installed on your system.\nPlease go to:\nhttps://computation.llnl.gov/casc/sundials/main.html\nand follow the instruction for installing SUNDIALS in MATLAB.');
else
%     disp('SUNDIALS Package FOUND!');
end

%Check for CasADi
% disp('Checking for CasADi availability...');

try
    import casadi.*
catch e
    error('The CasADi package does not seem to be present in the search path or installed on your system.\nPlease go to:\nhttps://github.com/casadi/casadi/wiki/InstallationInstructions\nand follow the instruction for installing CasADi in MATLAB.');
end

% disp('CasADi Package FOUND!');


% For each parameter structure, check if all the fields are set.
% for i=1:length(param)
%     %Check if the parameters have been set correctly
%     [result, missing] = checkBatteryParameters(param{i});
%     if(result==0)
%         disp('Warning, there are missing parameters in the param array.')
%         disp('Please find below, the list of missing parameters:')
%         for jj=1:length(missing)
%             disp(missing{jj});
%         end
%         error('Please fix the problem and restart the script')
%     end
% end

end
