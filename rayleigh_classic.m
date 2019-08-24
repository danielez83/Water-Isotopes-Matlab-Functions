% This function will compute the isotopic composition of Water Vapor from a
% water vapor mass of isotopic composition of d_0
% - d_0 starting composition of water vapor mass (delta units), e.g. -10‰
% - f fraction of water vapor remain
% - temperature (k), usually 270.15 K, mean adiabatically adjusted
% temperature at the lifted condensation level
% - isotope : 2 = D/H, 18 = 18O/16O

function rayleigh_classic = rayleigh_classic(d_0, f, temp, isotope)
    switch isotope
        case 2
            alfa = alpha2_LV(temp);
        case 18
            alfa = alpha18_LV(temp);
        otherwise
            disp('Bad Isotope value: 2 for D/H, 18 for 18O/16O')
    end
    rayleigh_classic = (((d_0/1000)+1)*f.^(alfa - 1) - 1)*1000; % Approximation from Kendal 1998, Fundamentals of isotope geochemistry
end
