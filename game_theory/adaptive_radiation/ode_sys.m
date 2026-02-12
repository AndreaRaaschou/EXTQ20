function [fy] = ode_sys( t , y , p )

%Competitors, GLV
% Set up the matrices of the GLV (see pp 348 in Case)
%See also Case pp 387 for an example of the full model with dispersal

% k=[p.r_A; p.mu_A];
k=p.g_d_vec;
A=p.A;
% disp_mat=p.disp_mat; 

N=y;

%fN=max(k.*N+(A*N).*N + disp_mat*N,0); %för att två consumenter med lågt sigmaK ska fungera 

%fN=k.*N+(A*N).*N + disp_mat*N; %för att consumenter och predatorerz ska fungera 

% fN=disp_mat*(k.*N+(A*N).*N);
fN=k.*N+(A*N).*N;

fy=fN;
