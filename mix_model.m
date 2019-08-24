% Isotopic composition of mixing water vapor
% This function compute the isotopic composition of a binary mixture from
% two end member charaterized by q_0 d_0 and q_1 d_1. 
% A complete description can be found in Noone, 2011.
% q = actual mixing ratio (mmol/mol, ppmv, ...)
% q_0 = mixing ratio of the lowest endmember (mmol/mol, ppmv, ...)
% d_0 = isotopic composition of the lowest end member (delta values)
% q_1 = mixing ratio of the higest endmember (mmol/mol, ppmv, ...)
% d_1 = isotopic composition of the higest end member (delta values)

function mix_model = mix_model(q,  q_0, d_0, q_1, d_1)
    dF = (q_1*d_1-q_0*d_0)/(q_1-q_0);
    mix_model = (q_0*(d_0-dF)*(1./q))+dF; % note the .* --> element by element
end
