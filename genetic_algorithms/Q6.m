%x = linspace(0, 10);
%y = linspace(0, 10);


%function [z] = sinus_fun(x, y)
%    z = x * sin(4 * x) + 1.1 * y * sin(2 * y);
%end


[X Y]=meshgrid(0:.1:10,0:.1:10);
Z=X.*sin(4*X)+1.1*Y.*sin(2*Y);
figure(1)
surf(X,Y,Z)

xlabel('x')
ylabel('y')
zlabel('z')
