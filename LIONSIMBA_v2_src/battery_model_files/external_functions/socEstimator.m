function Sout = socEstimator(t,t0,tf,states,extraData,param)
% SOCESTIMATOR estimates SOC of the cell.
% This function can be called after each integration step (within the time-stepping loop).
% It can be used to implement a custom SOC estimation by changing the function handle in
% Parameters_init suitably. It is mandatory to respect the sign of the function and to return
% a scalar value. The Sout value will be concatenated and it will be available in the result array
% at the end of the simulation.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%                   Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%                   Krishnakumar Gopalakrishnan and Ian Campbell, Imperial College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(param.SolidPhaseDiffusion~=3)
    cs_average = states.cs_average(end,param.Np+1:end);
else
    start_index = param.Nr_p*param.Np+1;
    end_index   = start_index+param.Nr_n-1;
    cs_average  = zeros(param.Nn,1);
    for n=1:param.Nn
        cs_average(n)   = 1/param.Rp_n*(param.Rp_n/param.Nr_n)*sum(states.cs_average(end,start_index:end_index));
        start_index     = end_index + 1;
        end_index       = end_index + param.Nr_n;
    end
end
% Estimates the SOC according to the current information of the states.
Csout= sum(cs_average);
Sout = 100*(1/param.len_n*(param.len_n/(param.Nn))*Csout/param.cs_maxn);

end
