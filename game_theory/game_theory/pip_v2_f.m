%% Invasion Fitness Landscape (MATLAB version)
% Pontarp exercise

clear; close all; clc;

%% Parameters
r   = 1;
mu0 = 0.2;
e   = 0.1;
a0  = 10;
k   = 5;
m   = 0.1;   % Try also m = 0

%% Trait values: resident x and mutant x'
x_vals = -1:0.005:1;

% Meshgrid for resident and mutant trait combinations
[Xres, Xmut] = meshgrid(x_vals, x_vals);

%% Attack rates
Ares = a0 - k .* (Xres.^2) / 2;
Amut = a0 - k .* (Xmut.^2) / 2;

%% Mortalities
MUres = mu0 + 0.5 .* m .* (Xres - 1).^2;
MUmut = mu0 + 0.5 .* m .* (Xmut - 1).^2;

%% Invasion fitness function
F = (MUres .* Amut ./ Ares) - MUmut;

%% Plot invasion fitness landscape
figure;
pcolor(Xres, Xmut, F);
shading interp;
colorbar;
caxis([-0.1 0.1]);

hold on;

% Contour where fitness = 0
contour(Xres, Xmut, F, [0 0], 'k', 'LineWidth', 2);

xlabel('Resident x');
ylabel("Invading x'");
title(['Invasion Fitness Landscape (m = ', num2str(m), ')']);

