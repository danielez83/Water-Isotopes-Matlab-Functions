% Fractionation factor of 2H Liquid/Vapor
% Ice-Liquid point at 268.15K (-5°C) 
% T [k]
function alpha2_LV = alpha2_LV(T)
    if T >= 268.15
        % Majoube 1971, Temperature in Kelvin
        alpha2_LV = exp((24844./(T.^2))-(76.248./T)+(0.052612));
    else
        alpha2_LV = exp((24844./(T.^2))-(76.248./T)+(0.052612)); % do not differentiate
        %alpha2_LV = exp((48888./(T.^2))-(203.10./T)+0.2133); % Ellehoj et al., 2013
    end
end
