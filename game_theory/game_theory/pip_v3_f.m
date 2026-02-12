%% Invasion Fitness Landscape (Two-resource consumer system)
% Translated from Python code (Pontarp)

clear; close all; clc;

%% Parameters
r  = 1;
mu = 0.3;
e  = 0.1;
a0 = 10;
k  = 1;

%% Trait grid
x_vals = -2:0.01:2;
[Xres, Xmut] = meshgrid(x_vals, x_vals);

%% Attack rates (resident and mutant on two resources)
A1res = a0 .* exp(-k .* (Xres + 1).^2);
A2res = a0 .* exp(-k .* (Xres - 1).^2);
A1mut = a0 .* exp(-k .* (Xmut + 1).^2);
A2mut = a0 .* exp(-k .* (Xmut - 1).^2);

%% Invasion fitness grid
F = zeros(size(Xres));

% Loop over resident-mutant combinations
for i = 1:size(Xres,1)
    for j = 1:size(Xres,2)

        a1r = A1res(i,j);
        a2r = A2res(i,j);
        a1m = A1mut(i,j);
        a2m = A2mut(i,j);

        % Solve resident equilibrium (R1, R2, C)
        [R1, R2, C, ok] = equilibrium(r, a1r, a2r, e, mu); %#ok<ASGLU>

        if ok
            % Mutant invasion fitness in resident environment
            F(i,j) = e * (a1m * R1 + a2m * R2) - mu;
        else
            F(i,j) = NaN; % handle singular / ill-conditioned case
        end
    end
end

%% Plotting
figure('Position',[200 200 800 600]);

% Filled contour plot (like contourf in Python)
contourf(Xres, Xmut, F, 50, 'LineColor', 'none');
colorbar;
ylabel(colorbar, 'Invasion Fitness');

hold on;
% Contour line where fitness = 0
contour(Xres, Xmut, F, [0 0], 'k', 'LineWidth', 1);

xlabel('Resident trait (x)');
ylabel("Mutant trait (x')");
title('Invasion Fitness Landscape');

axis tight;
box on;

%% --- Local function: equilibrium solver ---
function [R1, R2, C, ok] = equilibrium(r, a1, a2, e, mu)
% Solves the equilibrium resource-consumer system
%
% A * X = Y where X = [R1; R2; C]

    A = [r, 0, a1;
         0, r, a2;
         a1, a2, 0];

    Y = [r; r; mu/e];

    % Try to solve; handle singular / ill-conditioned matrices robustly
    ok = true;
    try
        % Equivalent to numpy.linalg.solve
        X = A \ Y;

        R1 = X(1);
        R2 = X(2);
        C  = X(3);

        % Optional sanity check (avoid crazy results if nearly singular)
        if any(~isfinite(X)) || rcond(A) < 1e-12
            ok = false;
            R1 = NaN; R2 = NaN; C = NaN;
        end

    catch
        ok = false;
        R1 = NaN; R2 = NaN; C = NaN;
    end
end
