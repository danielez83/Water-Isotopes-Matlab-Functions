% Approximation of Friction velocity from Wind Speed
% for open sea from Merlivat and Jouzel 1979
% result = ustar [cm/s]
% wind_speed [m/s]
% z, height [m]

function MJ79_ustar = MJ79_ustar(wind_speed, z)
    % http://www.energyhunters.it/content/misurare-il-vento-le-variazioni-del-vento-con-laltezza-dal-suolo
    z_0 = 3.5e-4; % roughness length (m)
    k = 0.4; % Von Karman constant, adimensional
    MJ79_ustar = 100*(wind_speed * k)/log(z/z_0); % [cm/s]
end
