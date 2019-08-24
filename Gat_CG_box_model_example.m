% This model compute the evolution of isotopic composition of water vapor 
% when an air parcel, with d_0 RH_0, is filled by local evaporation. 
% The starting air parcel could be also completely DRY 
% (d_0 = whatever you want because RH = 0).
% As reported in Gat et al., 2003 the RH and d_atmos are updated ech step.
% Se Gat at al., 2003 p958 for formulas. Here is used the full
% parametrization of Craig Gordon model and not the simplified.

% To run this script you need two matrixes as starting files:
% d18O_data --> [RH, Air_T, d18_Air, SST, d18_L] 
% dD_data --> [RH, Air_T, dD_Air, SST, dD_L] 
% example matrix are saved in "CG_box_model_example.mat"

% Load example data
load('CG_box_model_example.mat');

% Uncomment the following lines to change the isotopic composition of water [‰]
d18O_data(:,5) = 1.3%-2;
dD_data(:,5) = 6%-3.2;

% Uncomment the following lines to change air temperature [°C]
d18O_data(:,2) = 12;
dD_data(:,2) = d18O_data(:,2);

% Uncomment the following lines to change SST [°C]
d18O_data(:,4) = 12%10;
dD_data(:,4) = d18O_data(:,4);


%% Compute d18O
for t = 1:size(d18O_data,1)-1
    CG18 = CG_dE_18(d18O_data(t,1), (d18O_data(t,2)+273.15), d18O_data(t,3), (d18O_data(t,4)+273.15), d18O_data(t,5));
    d18O_data(t, 6) = CG18;
    d18O_data(t+1, 3) = ((d18O_data(t, 1)/100)*d18O_data(t, 3) + CG18*(d18O_data(t+1, 1) - d18O_data(t, 1))/100)/(d18O_data(t+1, 1)/100);
end

% Assign the same temperatures
dD_data(:,2) = d18O_data(:,2);
dD_data(:,4) = d18O_data(:,4);
%% Compute dD
for t = 1:size(dD_data,1)-1
    CG2 = CG_dE_2(dD_data(t,1), (dD_data(t,2)+273.15), dD_data(t,3), (dD_data(t,4)+273.15), dD_data(t,5));
    dD_data(t, 6) = CG2;
    dD_data(t+1, 3) = ((dD_data(t, 1)/100)*dD_data(t, 3) + CG2*(dD_data(t+1, 1) - dD_data(t, 1))/100)/(dD_data(t+1, 1)/100);
end

d_excess(1, 1) = dD_data(2,3) - 8 * d18O_data(2,3);
for t = 2:size(d18O_data,1)
    d_excess(t, 1) = dD_data(t,3) - 8 * d18O_data(t,3);
end

%% Show some nice plots
% Evolution of water vapor isotopic composition from RH_0 to RH=100%
figure
subplot(3,1,1)
plot(d18O_data(2:end, 1), d18O_data(2:end, 3))
ylabel('\delta^{18}O (‰)')
title('Evolution of \delta^{18}O')

subplot(3,1,2)
plot(dD_data(2:end, 1), dD_data(2:end, 3))
ylabel('\deltaD (‰)')
title('Evolution of \deltaD')

subplot(3,1,3)
plot(d18O_data(2:end, 1), d_excess(2:end, 1))
ylabel('d-excess (‰)')
title('Evolution of d-excess')

xlabel('RH (%)')



% Evolution on d18O vs dD plane with GMWL as a reference
figure
plot(d18O_data(2:end, 3), dD_data(2:end,3)) % Water vapor in the air parcel
xlabel('\delta^{18}O (‰)')
ylabel('\deltaD (‰)')
title('Evolution of WV isotopic composition')
% Compute GMWL
hold on 
axis_setup = gca;
GMWL = zeros(2,2);
GMWL(1,1) = axis_setup.XLim(1);
GMWL(2,1) = axis_setup.XLim(2);
GMWL(1,2) = axis_setup.XLim(1)*8+10;
GMWL(2,2) = axis_setup.XLim(2)*8+10;
plot(GMWL(:,1), GMWL(:,2), '--')
legend('WV air parcel','GMWL')
hold off
