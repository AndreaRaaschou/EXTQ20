% parameters:
r = 1; 
mu = 0.2;
e = 0.1;
a0 = 10;
k = 10;
 
% Create matrices with resident x (Xres) and mutant x (Xmut):
[Xres,Xmut] = meshgrid(-1:0.05:1);
% Calculate attack rates:
Ares = a0 - k*Xres.^2/2;
Amut = a0 - k*Xmut.^2/2;
 
% Invasion fitness:
F = mu*(Amut./Ares - 1);
 
% Plot:
figure
pcolor(Xres,Xmut,F)
caxis([-0.1 0.1])
shading interp
colorbar
hold on
contour(Xres,Xmut,F,[0 0],'k')
xlabel('Resident x')
ylabel('Invading x''')
