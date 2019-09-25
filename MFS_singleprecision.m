function [k_m,sigma_min]=MFS_singleprecision(colloc)
% Computes QR-MFS for ellipse with single precision

    format long;
% Parameters TO BE SET    
    % Initial guess for minimum
    init1=2.9;
    init2=init1;    %Only relevant in combination with fminbound instead of fminsearch, else set =init1
    % Interior points constant
    global m_I;
    m_I=10;
    % Number of sample functions
    global np;
    np=colloc;
    % Semi axis along x axis
    global a;
    a=1.0;    % Semi axis along y axis
    global b;
    b=1.0;
    % Interior points circle radius
    global R_i;
    R_i=0.5*min(a,b);
    % Exterior source points circle radius
    global R_e;
    R_e=5*max(a,b);
    % Index of refraction    
    global n;
    n=4;
    % Plot domain
    plot_flag=false;
    
%   Plot domain,normals, interior and exterior points
if (plot_flag==true)
    interval = linspace(0,2*pi);
    y1 = a*cos(interval);
    y2 = b*sin(interval);
    x1 = R_i*cos(interval);
    x2 = R_i*sin(interval);
    z1 = R_e*cos(interval);
    z2 = R_e*sin(interval);
    figure(1)
    plot(y1,y2);
    axis equal;
    hold on
    plot(x1,x2);
    hold on
    plot(z1,z2);
    hold on
    phi=linspace(0,2*pi,np+1);
    phi=phi(1:end-1); 
    quiver(a*cos(phi),b*sin(phi),b*cos(phi)./sqrt((b*cos(phi)).^2+(a*sin(phi)).^2),a*sin(phi)./sqrt((b*cos(phi)).^2+(a*sin(phi)).^2));
    title('ellipse & normal vectors');
    hold off
end   
 %   [k_m,sigma_min,exitflag]=fminbnd(@f,init1,init2,optimset('TolX',1e-16));
 %   [k_m,sigma_min,exitflag]=fminsearch(@f,init1,optimset('TolX',1e-16));
     fplot(@f,[0.1 3]);
end

function s=f(lambi)

    if(lambi==0)
        disp('Minimization approaches zero: Non-admissible wave number!');
        return
    end
    
    global m_I;
    global np
    global n;
    global a;
    global b;
    global R_e;
    global R_i;

    % create equidistant angles for scattering, artificial and interior parametrization
    phi_I=linspace(0,2*pi,m_I+1);
    phi_I=phi_I(1:end-1);
    interior_x=R_i*cos(phi_I);
    interior_y=R_i*sin(phi_I);
    interior=[interior_x;interior_y];    
    
    phi=linspace(0,2*pi,np+1);
    phi=phi(1:end-1); 
    
    sca_bou_x=a*cos(phi);
    sca_bou_y=b*sin(phi);
    scatter=[sca_bou_x;sca_bou_y];
    
    exterior_x=R_e*cos(phi);
    exterior_y=R_e*sin(phi);
    exterior=[exterior_x;exterior_y];
    
    normal_x=b*cos(phi)./sqrt((b*cos(phi)).^2+(a*sin(phi)).^2);
    normal_y=a*sin(phi)./sqrt((b*cos(phi)).^2+(a*sin(phi)).^2);
    normal=[normal_x;normal_y];

    % Setup up the interior Dirichlet boundary matrix 
    M_D_i=zeros(np,np);
    for i=1:np
      for j=1:np
        M_D_i(i,j)=kernelHelmholtz(scatter(:,i),exterior(:,j),lambi);
      end
    end      
    
    % Setup up the exterior Dirichlet boundary matrix 
    M_D_e=zeros(np,np);
    for i=1:np
      for j=1:np
        M_D_e(i,j)=kernelHelmholtz(scatter(:,i),exterior(:,j),sqrt(n)*lambi);
      end
    end  
    
    % Setup up the interior Neumann boundary matrix 
    M_N_i=zeros(np,np);
    for i=1:np
      for j=1:np
        M_N_i(i,j)=nabla_kernelHelmholtz(scatter(:,i),exterior(:,j),lambi)*normal(:,i);
      end
    end      
    
    % Setup up the exterior Neumann boundary matrix 
    M_N_e=zeros(np,np);
    for i=1:np
      for j=1:np
        M_N_e(i,j)=nabla_kernelHelmholtz(scatter(:,i),exterior(:,j),sqrt(n)*lambi)*normal(:,i);
      end
    end      
    
    % Setup up the interior points of interior domain 
    M_I_i=zeros(m_I,np);
    for i=1:m_I
      for j=1:np
        M_I_i(i,j)=kernelHelmholtz(interior(:,i),exterior(:,j),lambi);
      end
    end
    
    % Setup up the interior points of exterior domain 
    M_I_e=zeros(m_I,np);
    for i=1:m_I
      for j=1:np
        M_I_e(i,j)=kernelHelmholtz(interior(:,i),exterior(:,j),sqrt(n)*lambi);
      end
    end
    
    % Trefethen
    A=[M_D_i  M_D_e ; M_N_i  M_N_e ; M_I_i  zeros(m_I,np) ; zeros(m_I,np) M_I_e];
    %% COMMENT: ill-conditioning of A mostly affects QR, cf. condition number of QR factorization for (Trefethen's) accuracy estimate
    A=single(A);   
    [Q,~]=qr(A,0);
    tmp=svd(Q(1:(2*np),:)); 
    s=tmp(2*np);    % singular values stored in decreasing order
end

function z=kernelHelmholtz(x,y,k)
    r=norm(x-y,2);
    z=besselh(0,1,k*r);
end

function z=nabla_kernelHelmholtz(x,y,k) 
    r=norm(x-y,2);
    z=-k*besselh(1,1,k*r)*(x-y)'/r;
end
