function  [k_m,sigma_min]=MFS_ITE_ellipse_anisotropic(m,init,semi_x,semi_y,Rrel_S,n,a_xx,a_yy,a_xy)
% ROUTINE: Computes approximate ITE for ellipse with modified MFS
% INPUT:  
%   m = number of boundary collocation points (or equivalently source points)
%   init = initial guess for approximate ITE location
%   semi_x = semi axis of ellipse in x direction
%   semi_y = semi axis of ellipse in y direction    
%   Rrel_S = scaling factor of radius for circular source point radius wrt semi major axis
%   n = refractive index
%   a_xx = a11 of symmetric anisotropy matrix
%   a_yy = a22 of symmetric anisotropy matrix
%   a_xy = a12 of symmetric anisotropy matrix

    format long
    
% Interior points parameters
    % Number of interior points
    m_I=10;
    % Scaling factor wrt minor axis for radius R_I of circular interior points
    Rrel_I=0.5;
    % Absolute radius for inner and source points
    R_I=Rrel_I*min(semi_x,semi_y);
    R_S=Rrel_S*max(semi_x,semi_y);
    
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
    f_tmp=@(kappa) (f(kappa,m,m_I,interior,boundary,sources,normals,n,a_xx,a_yy,a_xy)); 
    [k_m,sigma_min]=fminsearch(f_tmp,init,optimset('TolX',1e-16));
    %[k_m,sigma_min]=fminbnd(f_tmp,init,init+0.5,optimset('TolX',1e-16));
          
%%   Plot scattering domain, boundary normals, interior and source points
%     figure(m)
%     fplot(f_tmp, [3.338 3.339]);
%
%     interval = linspace(0,2*pi);
%     y1 = semi_x*cos(interval);
%     y2 = semi_y*sin(interval);
%     x1 = R_I*cos(interval);
%     x2 = R_I*sin(interval);
%     z1 = Rrel_S*y1;
%     z2 = Rrel_S*y2;
%     figure(2)
%     plot(y1,y2);
%     axis equal;
%     hold on
%     plot(x1,x2);
%     hold on
%     plot(z1,z2);
%     hold on
%     phi=linspace(0,2*pi,m+1);
%     phi=phi(1:end-1); 
%     quiver(semi_x*cos(phi),semi_y*sin(phi),semi_y*cos(phi)./sqrt((semi_y*cos(phi)).^2+(semi_x*sin(phi)).^2),semi_x*sin(phi)./sqrt((semi_y*cos(phi)).^2+(semi_x*sin(phi)).^2));
%     hold off
%     title('ellipse');
end

function s=f(k,m,m_I,pts_I,pts_B,pts_S,normals,n,a_xx,a_yy,a_xy)

    if(k==0)
        s=0;
        return
    end
    
    % Unit vectors for partial derivatives
    e_x=[1 ; 0];
    e_y=[0 ; 1];
    
    % Anisotropic coefficients
    A=[a_xx a_xy ; a_xy a_yy];
    A_root = sqrtm(A);
 
    % Setup up Dirichlet boundary for v
    D_v=zeros(m,m);
    for i=1:m
      for j=1:m
        D_v(i,j)=kernelHelmholtz(pts_B(:,i),pts_S(:,j),k);
      end
    end      
    
    % Setup up Dirichlet boundary for w
    D_w=zeros(m,m);
    for i=1:m
      for j=1:m
        D_w(i,j)=kernelHelmholtz(A_root\pts_B(:,i),A_root\pts_S(:,j),sqrt(n)*k);
      end
    end  
    
    % Setup up Neumann boundary for v
    N_v=zeros(m,m);
    for i=1:m
      for j=1:m
        N_v(i,j)=nabla_kernelHelmholtz(pts_B(:,i),pts_S(:,j),k)*normals(:,i);
      end
    end      
    
    % Setup up Dirichlet boundary for w
    N_w=zeros(m,m);
    for i=1:m
      for j=1:m
        N_w(i,j)=nabla_kernelHelmholtz(A_root\pts_B(:,i),A_root\pts_S(:,j),sqrt(n)*k)*A_root*normals(:,i);
      end
    end      
    
    % Setup up interior samples of v
    I_v=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        I_v(i,j)=kernelHelmholtz(pts_I(:,i),pts_S(:,j),k);
      end
    end
    % Setup up partial-x interior samples of v 
    I_vx=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        I_vx(i,j)=nabla_kernelHelmholtz(pts_I(:,i),pts_S(:,j),k)*e_x;
      end
    end
    % Setup up partial-y interior samples of v 
    I_vy=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        I_vy(i,j)=nabla_kernelHelmholtz(pts_I(:,i),pts_S(:,j),k)*e_y;
      end
    end
    
    % Setup up interior samples of w
    I_w=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        I_w(i,j)=kernelHelmholtz(A_root\pts_I(:,i),A_root\pts_S(:,j),sqrt(n)*k);
      end
    end
    % Setup up partial-x interior samples of w 
    I_wx=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        I_wx(i,j)=nabla_kernelHelmholtz(A_root\pts_I(:,i),A_root\pts_S(:,j),sqrt(n)*k)*A_root*e_x;
      end
    end
    % Setup up partial-y interior samples of v
    I_wy=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        I_wy(i,j)=nabla_kernelHelmholtz(A_root\pts_I(:,i),A_root\pts_S(:,j),sqrt(n)*k)*A_root*e_y;
      end
    end
    
    % Trefethen QR utilities
    A=[D_v  D_w ; N_v  N_w ; I_v  zeros(m_I,m) ; I_vx  zeros(m_I,m); I_vy  zeros(m_I,m); zeros(m_I,m) I_w; zeros(m_I,m) I_wx; zeros(m_I,m) I_wy];
    [Q,~]=qr(A,0);
    tmp=svd(Q(1:(2*m),:));
    s=tmp(2*m);    
end

function z=kernelHelmholtz(x,y,k)
    r=norm(x-y,2);
    z=besselh(0,1,k*r);
end

function z=nabla_kernelHelmholtz(x,y,k) 
    r=norm(x-y,2);
    z=-k*besselh(1,1,k*r)*(x-y)'/r;
end
