% Estimate isotopic composition of water vapor flux with
% the Craig-Gordon model (for 18O)
% Computation is made for sea water (activity 0.98), see Horita 2008 for
% additional informations
% RH: 1 - 100 [%]
% Ait temperature [k]
% d18O of atmosphere [‰] 
% Sea/Water temperature [k]
% d18O of water [‰]
% For the computation the following assumption are made:
% - RH is not normalized to SST, it will be normalized to,
% - Theta parameter = 1 (Gat 1996)
% - n = 1/2, exponent for diffusion parameter is for ocean (Gat 1996)
% - kinetic fract C_k as defined in Merlivat (1978)

%------------------------------------------------------------------------------------
% Notes on saturation water vapor pressure. The formula used in this
% algorithm is based on the Engineering Toolbox one
% (http://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html).
% Reference for this formula can be taken from:
% - Vladilo et al. 2013, "The habitable zone of Earth-like planets with different
%   levels of atmospheric pressure"
% - Rao et al. 2008, "Convective condensation of vapor in the presence of a
%   non-condensable gas of high concentration in laminar flow in a vertical pipe"
% - Saraireh & Thorpe 2010, "Modelling of Heat and Mass Transfer Involvin
%   vapour ..."
%------------------------------------------------------------------------------------

% % HORITA VERSION
function CG_dE_18 = CG_dE_18(RH, T_air, d18_atmos, T_water, d18_water)

    % For THETA chose one of the following values:
    % 0.5 Mediterranean Sea
    % 0.88 great lake
    % 1 if you don't know...
    % THETA = 0.67; % Ok for Venice
     THETA = 0.52; % continuous measurements/open sea
    % THETA = 1; % for soil evaporation
    
    % For activity chose one of the following values:
    % 0.98 Sea water
    % 1 Freshwater
    activity = 0.98;
    % activity = 1;
    
    % For n chose one of the following values:
    % 1/2 for ocean, 0.5 fully tubolent conditions or rough surfaces,
    % 2/3 smooth surfaces (Brusaert 1982)
    % 1 stagnant layer (Gat 1996)
    % n = 0.48; % Ok for Venice
    n = 0.5;
    %n = 1; % for soil evaporation
    
    % normalize humidity at the surface 
    e_a = (RH./100)*(exp(77.3450 + 0.0057 .* (T_air) - (7235./(T_air)))./(T_air).^8.2);
    e_s = (exp(77.3450 + 0.0057 .* (T_water) - (7235./(T_water)))./(T_water).^8.2);
    one_minus_h = (activity*e_s-e_a)/(activity*e_s);
    
    C_k = n * (1 - 0.9727) * 1000; % 0.9691 (Cappa et al. 2003); 0.9727 (Merlivat 1978)
    alpha18_VL = 1./alpha18_LV(T_water);
    % epsilon_k = THETA * (1 - (activity*norm_RH(RH, T_air, T_water)/100)) * C_k;
    epsilon_k = THETA * one_minus_h * C_k;
    epsilon_star = (1 - alpha18_VL)*1000;
    
    CG_dE_18 = ((alpha18_VL * d18_water) - ((1-one_minus_h) * d18_atmos) - (epsilon_star + epsilon_k))/(one_minus_h + (0.001 * epsilon_k));
end
