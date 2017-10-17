function [resDfilm, rhsDfilm] = SEI_layer(dfilmThickness, jflux_side, param)
% SEI_LAYER  Describes the dynamics of the Solid-Electrolyte layer at the anode side.

% This code was written by Marcello Torchio, University of Pavia.
% Please send comments or questions to
% marcello.torchio01@ateneopv.it
%
% Copyright 2015: 	Marcello Torchio, Lalo Magni, and Davide M. Raimondo, University of Pavia
%					Bhushan Gopaluni, University of British Columbia
%                 	Richard D. Braatz, MIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhsDfilm = ( -jflux_side * param.M_n / (param.rho_n));
resDfilm = dfilmThickness - rhsDfilm;

end