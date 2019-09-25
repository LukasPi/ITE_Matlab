function exact_ITE_unitdisc_isotropic(order)
% Computes 2D ITEs of unit disc and their radial-symmetric eigenfunctions from Fourier-Bessel ansatz of corresponding order
format long
% Index of refraction
n=4;%+0.5*i;
% Neighbourhood of ITE
x_0=13; 

% Use: J'(m)=0.5*(J(m-1)-J(m+1)
Dbesselj=@(x)(0.5*(besselj(order-1,x)-besselj(order+1,x)));
% Determinant normalized by trivial root to avoid flattening of graph
det=@(z) (10^order*(besselj(order,z)*sqrt(n).*Dbesselj(sqrt(n)*z)-besselj(order,sqrt(n)*z).*Dbesselj(z))./(z/2).^(2*order-1));
%Plot determinant function 
figure(1)
    fplot(det, [x_0-1 x_0+1])
    hold off
figure(2)
    [X,Y] = meshgrid(0:0.01:4,-1:0.01:1);
    H =abs(det(X+1i*Y));
    contour(X,Y,H,0:0.01:0.1)
    hold off
 
 %% Compute exact ITE via vanishing determinant
 options = optimoptions('fsolve','TolFun',1e-16);
 k=fsolve(det,x_0,options) %aborts if det(x_0) is relatively small, but no root
 %k=fzero(det,x_0) %real-valued k

%% Compute L2-norms of eigenfunctions
% Bessel coefficients c for v in equation (10) of paper: "c*v-w"=0 at r=1
c=besselj(order,sqrt(n)*k)/besselj(order,k);
% Compute L2-integral of v and w by Fubini wrt radial(r) and angular(a) part (note 2D Polar-Jacobian is r):
    % ...square as absolute values
        %v_r = @(r) (abs(c^2*besselj(order,k*r).^2.*r*pi));
        %w_r = @(r) (abs(besselj(order,sqrt(n)*k*r).^2.*r*pi));
    % ...complex square without positivity
        v_r = @(r) (c^2*besselj(order,k*r).^2.*r*pi);
        w_r = @(r) (besselj(order,sqrt(n)*k*r).^2.*r*pi);
   
 L2_v=integral(v_r,0,1,'AbsTol',1e-15);
 L2_w=integral(w_r,0,1,'AbsTol',1e-15);
 temp=L2_v+L2_w;
 L2_v=L2_v/temp;
 L2_w=L2_w/temp;
 critical=(L2_v-n*L2_w)
 
 end