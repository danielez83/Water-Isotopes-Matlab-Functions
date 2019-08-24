% Estimate isotopic composition of water vapor flux with
% the Craig-Gordon model (for 2H) MJ79
% Computation is made for sea water (activity 0.98), see Horita 2008 for
% additional informations
% RH: 1 - 100 [%]
% Ait temperature [k]
% dD of atmosphere [‰]
% Sea/Water temperature [k]
% dD of water [‰]
% Wind dep 'on' or 'off', if off 3.9‰ will be used folowing Pfahl & Wernli (2009)
% Wind_Speed [m/s] as described in Merlivat and Jouzel (1979)
% For the computation the following assumption are made:
% - RH is not normalized to SST, it will be normalized to,

%------------------------------------------------------------------------------------
% Notes on saturation water vapor pressure. The formula used in this
% algorithm is based on the one found on Engineering Toolbox
% (http://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html).
% Reference for this formula can be taken from:
% - Vladilo et al. 2013, "The habitable zone of Earth-like planets with different
%   levels of atmospheric pressure"
% - Rao et al. 2008, "Convective condensation of vapor in the presence of a
%   non-condensable gas of high concentration in laminar flow in a vertical pipe"
% - Saraireh & Thorpe 2010, "Modelling of Heat and Mass Transfer Involvin
%   vapour ..."
%------------------------------------------------------------------------------------

function CG_dE_2_MJ79 = CG_dE_2_MJ79(RH, T_air, d2_atmos, T_water, d2_water, wind_dep ,Wind_Speed)
    % Constants
    R2SMOW = 155.76e-6;    % VSMOWS 2H Absolute concentration
    k = 7.5e-3*0.88;         % In case of no wind dependance
    activity = 0.98;        %Activity of the water, Freshwater = 1, Seawater = 0.98
    Wind_Speed_Threshold = 7; %m/s @ 10m height
    z_measure = 600; % wind speed height measure in cm
       
    % Normalize humidity
    e_a = (RH/100)*(exp(77.3450 + 0.0057 * (T_air) - (7235/(T_air)))/(T_air)^8.2);
    e_s = (exp(77.3450 + 0.0057 * (T_water) - (7235/(T_water)))/(T_water)^8.2);
    one_minus_h = (activity*e_s-e_a)/(activity*e_s);
    RH0 = (1 - one_minus_h)*100;
       
    % Estimate k-value
    if strcmp(wind_dep, 'on')
        if Wind_Speed > Wind_Speed_Threshold % rough regime
            ustar = MJ79_ustar(Wind_Speed, 10);
            k = 0.001*MJ79_k(ustar, z_measure, 2, 'rough');
        else % smooth regime
            ustar = MJ79_ustar(Wind_Speed, 10);
            k = 0.001*MJ79_k(ustar, z_measure, 2, 'smooth');
        end
    end
    
    % USE R version of CRAIG GORDON MODEL instead of delta version -------------
    h = RH0/100;
    Rl = ((d2_water/1000)+1)*R2SMOW;
    Ra = ((d2_atmos/1000)+1)*R2SMOW;
    Re = (1-k)*(((1/alpha2_LV(T_water))*Rl)-(h*Ra))/(1-h);
    CG_dE_2_MJ79 = ((Re/R2SMOW)-1)*1000;
end
