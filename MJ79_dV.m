% Isotopic composition of water vapor with the the MJ79 closure assumption
% result = [d18O, dD]
% SST [k]
% RH_SST [%]
% Wind Regime 3 or 8 [m/s]  
% dL = isotopic composition of water [‰]

function MJ79_dV = MJ79_dV(SST, RH_SST, regime, dL_18, dL_2)
    switch regime
        case 3
            k18 = 6.4;
        case 8
            k18 = 3.5;
    end
    d18O = (((1/alpha18_LV(SST))*(dL_18+1)*(k18))/(1-(k18*RH_SST/100)))-1;

    dD = (((1/alpha2_LV(SST))*(dL_2+1)*((k18*0.88)))/(1-(k18*0.88*RH_SST/100)))-1;
    MJ79_dV = [d18O, dD];
end
