function [k_m,sigma_min]=test_MFSstandard(m,init,semi_x,semi_y,Rrel_S,n)
% ROUTINE: Computes approximate ITE for ellipse with standard MFS
% INPUT:  
%   m = number of boundary collocation points (or equivalently source points)
%   init = initial guess for approximate ITE location
%   semi_x = semi axis of ellipse in x direction
%   semi_y = semi axis of ellipse in y direction    
%   Rrel_S = scaling factor of radius for circular source point radius wrt semi major axis
%   n = refractive index

% Parameters TO BE SET    
    % Interior points parameters
    % Number of interior points
    m_I=10;
    % Scaling factor wrt minor axis for radius R_I of circular interior points
    Rrel_I=0.5;
    % Absolute radius for inner and source points
    R_I=Rrel_I*min(semi_x,semi_y);
    R_S=Rrel_S*max(semi_x,semi_y);
    % Plot domain
    plot_flag=false;
    
    format long;

% Generation of computational points 
    % Interior points 
    phi_I=linspace(0,2*pi,m_I+1);
    phi_I=phi_I(1:end-1);
    interior_x=cos(phi_I);
    interior_y=sin(phi_I);
    interior=[R_I*interior_x ; R_I*interior_y]; 
    
    % Boundary collocation points
    phi=linspace(0,2*pi,m+1);
    phi=phi(1:end-1); 
    boundary=[semi_x*cos(phi); semi_y*sin(phi)];
    
    % Normal vectors along boundary
    normal_x=semi_y*cos(phi)./sqrt(semi_y^2*cos(phi).^2+semi_x^2*sin(phi).^2);
    normal_y=semi_x*sin(phi)./sqrt(semi_y^2*cos(phi).^2+semi_x^2*sin(phi).^2);
    normals=[normal_x ; normal_y];
    
    % Source points 
    sources=[R_S*cos(phi) ; R_S*sin(phi)];

% Modified MFS main routine
    f_tmp=@(kappa) (f(kappa,m,m_I,interior,boundary,sources,normals,n)); 
    %[k_m,sigma_min]=fminsearch(f_tmp,init,optimset('TolX',1e-16));
    %[k_m,sigma_min]=fminbnd(f_tmp,init,init+0.5,optimset('TolX',1e-16));
    
     figure(m)
     fplot(f_tmp, [1.8 4.7]);

%   Plot domain,normals, interior and exterior points
    if (plot_flag==true)
        interval = linspace(0,2*pi);
        y1 = semi_x*cos(interval);
        y2 = semi_y*sin(interval);
        x1 = R_I*cos(interval);
        x2 = R_I*sin(interval);
        z1 = R_S*cos(interval);
        z2 = R_S*sin(interval);
        figure(1)
        plot(y1,y2);
        axis equal;
        hold on
        plot(x1,x2);
        hold on
        plot(z1,z2);
        hold on
        phi=linspace(0,2*pi,m+1);
        phi=phi(1:end-1); 
        quiver(semi_x*cos(phi),semi_y*sin(phi),semi_y*cos(phi)./sqrt((semi_y*cos(phi)).^2+(semi_x*sin(phi)).^2),semi_x*sin(phi)./sqrt((semi_y*cos(phi)).^2+(semi_x*sin(phi)).^2));
        title('ellipse & normal vectors');
        hold off
    end   
end

function s=f(lambi,m,m_I,pts_I,pts_B,pts_S,normals,n)

    if(lambi==0)
        s=0;
        return
    end

    % Setup up the interior Dirichlet boundary matrix 
    M_D_i=zeros(m,m);
    for i=1:m
      for j=1:m
        M_D_i(i,j)=kernelHelmholtz(pts_B(:,i),pts_S(:,j),lambi);
      end
    end      
    
    % Setup up the exterior Dirichlet boundary matrix 
    M_D_e=zeros(m,m);
    for i=1:m
      for j=1:m
        M_D_e(i,j)=kernelHelmholtz(pts_B(:,i),pts_S(:,j),sqrt(n)*lambi);
      end
    end  
    
    % Setup up the interior Neumann boundary matrix 
    M_N_i=zeros(m,m);
    for i=1:m
      for j=1:m
        M_N_i(i,j)=nabla_kernelHelmholtz(pts_B(:,i),pts_S(:,j),lambi)*normals(:,i);
      end
    end      
    
    % Setup up the exterior Neumann boundary matrix 
    M_N_e=zeros(m,m);
    for i=1:m
      for j=1:m
        M_N_e(i,j)=nabla_kernelHelmholtz(pts_B(:,i),pts_S(:,j),sqrt(n)*lambi)*normals(:,i);
      end
    end      
    
    
    % Trefethen
    A=[M_D_i  M_D_e ; M_N_i  M_N_e];
    tmp=svd(A); 
    s=tmp(2*m);    % singular values stored in decreasing order
end

function z=kernelHelmholtz(x,y,k)
    r=norm(x-y,2);
    z=besselh(0,1,k*r);
end

function z=nabla_kernelHelmholtz(x,y,k) 
    r=norm(x-y,2);
    z=-k*besselh(1,1,k*r)*(x-y)'/r;
end
