function z=exact_ITE_disc_elastic
% Computes first exact 2D ITEs for unit disc from Sun paper (NOT all though since eigenfunctions are not always rotational symmetric!!!)
% Only besselj(1,:) feasible as required to solve Bessel equation of order 1

format long;

% Lamé constants
u=1/16;
l=1/4;
% densities
rho=4;
% radius of disc
R=0.5;
% Neighbourhood of ITE
x_0=1.45;  

% Abbreviations
a_1=sqrt(1/u);%sqrt(1/(2*u+l));
a_2=sqrt(rho/u);%sqrt(rho/(2*u+l));

% Use: J'(m)=0.5*(J(m-1)-J(m+1)
Dbesselj=@(x)(0.5*(besselj(0,x)-besselj(2,x)));
% Find roots of determinant function: 
det=@(z) (a_2*besselj(1,a_1*R*z).*Dbesselj(a_2*R*z)-a_1*besselj(1,a_2*R*z).*Dbesselj(a_1*R*z));
options = optimoptions('fsolve','TolFun',1e-17);
k=fsolve(det,x_0,options) %k=fzero(det,x_0) for real-valued k
%det(k)

[X,Y] = meshgrid(0:0.01:3,-2:0.01:2);
 H =det(X+1i*Y);
 contour(X,Y,abs(H),0:0.1:1)
end