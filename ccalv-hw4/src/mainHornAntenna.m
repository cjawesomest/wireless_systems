% main.m
%     Cameron J. Calv
%     
% Demo.m
%     Kapil R. Dandekar
%     Xaime Rivas Rey
%	  ECE-T512 - Wireless Systems Matlab Courseware
%        
%     This program illustrates the radiation intensity of a horned antenna
%% 0. The antenna parameters
clear all; clc; close all;

min_theta = 1*10^(-5);
max_theta = 2*pi-1*10^(-5);
theta_steps = 200;
thetas = linspace(min_theta, max_theta, 100);

min_phi = -pi/2;
max_phi = pi/2;
phi_steps = 200;
phis = linspace(max_phi, min_phi, phi_steps);

frequency = 2.4*10^9;
lambda = (3*10^8)/frequency; 
E_zero = 1;
r = 1;
% r = 1*10^(-7);
k = (2*pi)/lambda;
A = 0.5*lambda;
B = 0.25*lambda;
a = 0.25*lambda;
b = 0.25*lambda;
eta_not = (1.26*10^(-6))*(3*10^8);
efficiency = 0.8;

%% 1. Calculations in 3D

Zin = [];
integral = 0;
for phi = phis
    theta_row = [];
    for theta = thetas
        propagation_x = (A/lambda)*sin(theta)*cos(phi);
        propagation_y = (B/lambda)*sin(theta)*sin(phi);
        field_initial_1 = (4/pi)*((cos(pi*propagation_x))/...
            (1-4*(propagation_x^2)));
        field_initial_0 = (2/pi)*((sin(pi*propagation_y))/...
            (propagation_y));
        E_field_theta = j*(exp(-1*j*k*r)/(lambda*r))*E_zero*((A*B)/4)*...
            ((1+cos(theta))/2)*sin(phi)*field_initial_1*field_initial_0;
        E_field_phi = j*(exp(-1*j*k*r)/(lambda*r))*E_zero*((A*B)/4)*...
            ((1+cos(theta))/2)*cos(phi)*field_initial_1*field_initial_0;
        E_field_theta_mag = sqrt(real(E_field_theta)^2+imag(E_field_theta)^2);
        E_field_phi_mag = sqrt(real(E_field_phi)^2+imag(E_field_phi)^2);
        radiation_intensity = (r^2)*(E_field_theta_mag^2 + E_field_phi_mag^2)/(2*eta_not);
        theta_row = [theta_row, 10*log10(radiation_intensity)];
        integral = integral + radiation_intensity*sin(theta);
    end
    Zin = [Zin; theta_row];
end

intensities_mW = 10.^(Zin./10);
directivity = max(max(intensities_mW))/(integral);
gain = 10*log10(efficiency*directivity); %130dB

%% 2. Plots!
close all;
figure(1)
[Zs, Xs, Ys] = sphere3d(Zin,min_theta,max_theta,min_phi,max_phi, r, 1,'surf');
title("Horned Antenna Radiation Pattern");
xlabel("Z", 'Visible', 'on');
ylabel("X", 'Visible', 'on');
zlabel("Y", 'Visible', 'on');
% set(gca, 'visible', 'on')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
box on

azimuth_angle = 3*pi/2; %Theta
theta_idx = find((thetas>=azimuth_angle) == 1,  1);
azimuth_Xs = Xs(:, theta_idx)';
azimuth_Ys = Ys(:, theta_idx)';
figure(2);
plot(azimuth_Xs, azimuth_Ys);
set(gca, 'Visible', 'off')
title("Azimuth Plane", 'Visible', 'on');
xlabel("X", 'Visible', 'on');
ylabel("Y", 'Visible', 'on');

elevation_angle = pi/2; %Phi
phi_idx = find((phis>=elevation_angle) == 1,  1);
elevation_Ys = Xs(phi_idx, :);
elevation_Zs = -1*Zs(phi_idx, :);
figure(3);
plot(elevation_Ys, elevation_Zs);
set(gca, 'Visible', 'off')
title("Elevation Plane", 'Visible', 'on');
xlim([-2.1*10^(-17), 2.1*10^(-17)]);
xlabel("Y", 'Visible', 'on');
ylabel("Z", 'Visible', 'on');



