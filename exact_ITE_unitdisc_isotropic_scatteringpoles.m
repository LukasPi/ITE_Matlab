function exact_ITE_unitdisc_isotropic_scatteringpoles(order)
% Computes scattering poles for unit disc from Fourier-Hankel ansatz (~radiating) of corresponding order

format long
% Index of refraction
n=4;
% Neighbourhood of ITE
x_0=0.5-1.8i; 

% Use: H'(m)=0.5*(H(m-1)-H(m+1)
Dbesselh=@(x)(0.5*(besselh(order-1,x)-besselh(order+1,x)));
det=@(z) (besselh(order,z)*sqrt(n).*Dbesselh(sqrt(n)*z)-besselh(order,sqrt(n)*z).*Dbesselh(z));
% Plot determinant function 
 [X,Y] = meshgrid(0:0.01:5,-5:0.01:5);
 H =det(X+1i*Y);
 contour(X,Y,abs(H),0:0.01:1)

 %% Compute exact ITE via vanishing determinant
 %options = optimoptions('fsolve','TolFun',1e-16);
 %k=fsolve(det,x_0,options)
 %det(k)
end