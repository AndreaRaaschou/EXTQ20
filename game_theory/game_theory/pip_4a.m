% 4a
% Plot the Pairwise Invasibility Plot (PIP) of the
% Lotka-Volterra predator-prey model with logistic growth

clear; close all; clc;

r = 1;       % intrinsic growth rate
mu = 0.2;    % consumer mortality
eps = 0.1;   % conversion coefficient
a0 = 10;     % maximum attack rate (when x = 0)
k = 10;      % constans?

% use meshgrid to create a 2D matrix with (columns, rows)
x = linspace(-1, 1, 200);
y = linspace(-1, 1, 200); 
[X, Y] = meshgrid(x, y);

% calculate fitness for invading x' (Y)
Z = mu*(a(Y, a0, k)./a(X, a0, k) - 1);

% plot
figure
pcolor(X, Y,Z)
shading interp
colorbar
hold on
contour(X, Y, Z, [0, 0], 'k')
xlabel('Resident x', 'FontSize', 16)
ylabel('Invading x''', 'FontSize', 16)
hold off

% calculate attack rate
function [attack_rate] = a(x, a0, k)
attack_rate = a0 - 1/2 * k*x.^2;
end