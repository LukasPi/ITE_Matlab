function  [k_m,sigma_min]=MFS_ITE_Lshape_anisotropic(m,init,Rrel_S,n,a_xx,a_yy,a_xy)    %MFS_ITE_Lshape_anisotropic(105,10.1,1.01,1,1/4,1/4,0)
% ROUTINE: Computes approximate ITE for L-shape of unit length with modified MFS
% INPUT:  
%   m = number of boundary collocation points (or equivalently source points)
%   init = initial guess for approximate ITE location    
%   Rrel_S = scaling factor for scaling-translated sources wrt boundary
%   n = refractive index
%   a_xx = a11 of symmetric anisotropy matrix
%   a_yy = a22 of symmetric anisotropy matrix
%   a_xy = a12 of symmetric anisotropy matrix

    format long;

% Interior points parameters
    p=4;
    % Number of interior points
    m_I=20;
    % Scaling factor wrt inner circumference 
    Rrel_I=0.3;
    % Absolute radius for inner points (in 3rd quadrant)
    R_I=Rrel_I*0.5/tan(2*pi/(2*p));
    
% Generation of computational points 
    % Interior points 
    phi_I=linspace(0,2*pi,m_I+1);
    phi_I=phi_I(1:end-1);
    interior_x=cos(phi_I);
    interior_y=sin(phi_I);
    interior=[R_I*interior_x-0.25 ; R_I*interior_y-0.25]; 
    
    % Boundary collocation points
        % midpoint to edge
        ME=0.5/tan(pi/p);
        % midpoint to vertex
        MV=0.5/sin(pi/p);
        % generate boundary coordinates
        boundary_x=zeros(1,m);
        boundary_y=zeros(1,m);
        normal_x=zeros(1,m);
        normal_y=zeros(1,m);
        boundary_x(1,1)=ME;
        boundary_y(1,1)=0;
        normal_x(1)=1.;
        normal_y(1)=0.;
        for ii=1:(m-1)
            section=floor(ii*p/m+0.5);
            alpha=2*pi*ii/m-((section-1)*2*pi/p+pi/p);
            R=MV*cos(-pi/p)/cos(-pi/p+alpha);
            if(ii/m<=0.25)
                boundary_x(1,ii+1)=0.5-R*sin(2*pi*ii/m);
                boundary_y(1,ii+1)=0.5-R*cos(2*pi*ii/m);
                normal_x(1,ii+1)=sin(section*2*pi/p);
                normal_y(1,ii+1)=cos(section*2*pi/p);
            else  
                boundary_x(1,ii+1)=R*cos(2*pi*ii/m );
                boundary_y(1,ii+1)=R*sin(2*pi*ii/m );
                normal_x(1,ii+1)=cos(section*2*pi/p);
                normal_y(1,ii+1)=sin(section*2*pi/p);
            end
        end
    boundary=[boundary_x ; boundary_y];
%      phi=linspace(0,2*pi,m+1);
%      phi=phi(1:end-1);
%      sources=[Rrel_S*MV*cos(phi) ; Rrel_S*MV*sin(phi)];
    sources=[Rrel_S*boundary_x+0.25*(Rrel_S-1) ; Rrel_S*boundary_y+0.25*(Rrel_S-1)];
    normals=[normal_x ; normal_y];
    
    % Modified MFS main routine
    f_tmp=@(kappa) (f(kappa,m,m_I,interior,boundary,sources,normals,n,a_xx,a_yy,a_xy)); 
    %[k_m,sigma_min]=fminsearch(f_tmp,init,optimset('TolX',1e-16));
    %[k_m,sigma_min]=fminbnd(f_tmp,init-0.1,init+0.1,optimset('TolX',1e-16));
    
%%   Plot scattering domain, boundary normals, interior and source points  
      figure(m)
      fplot(f_tmp, [0 8])
      title('Lshape')

    % midpoint to edge
    ME=0.5/tan(pi/p);
    % midpoint to vertex
    MV=0.5/sin(pi/p);
    % generate boundary coordinates and their normal
    boundary_x=zeros(m ,1);
    boundary_y=zeros(m ,1);
    normal_x=zeros(m ,1);
    normal_y=zeros(m ,1);
    boundary_x(1)=ME;
    boundary_y(1)=0;
    normal_x(1)=1.;
    normal_y(1)=0.;
    figure(2)
    quiver(boundary_x(1),boundary_y(1),normal_x(1),normal_y(1));
    axis equal;
    hold on
    for ii=1:(m -1)
        section=floor(ii*p/m +0.5);
        alpha=2*pi*ii/m -((section-1)*2*pi/p+pi/p);
        R=MV*cos(-pi/p)/cos(-pi/p+alpha);
            if(ii/m<=0.25)
                boundary_x(ii+1)=0.5-R*sin(2*pi*ii/m);
                boundary_y(ii+1)=0.5-R*cos(2*pi*ii/m);
                normal_x(ii+1)=sin(section*2*pi/p);
                normal_y(ii+1)=cos(section*2*pi/p);
            else  
                boundary_x(ii+1)=R*cos(2*pi*ii/m );
                boundary_y(ii+1)=R*sin(2*pi*ii/m );
                normal_x(ii+1)=cos(section*2*pi/p);
                normal_y(ii+1)=sin(section*2*pi/p);
            end
        hold on
        quiver(boundary_x(ii+1),boundary_y(ii+1),normal_x(ii+1),normal_y(ii+1));
        hold on
    end
    interval = linspace(0,2*pi);
    x1 = R_I*cos(interval)-0.25;
    x2 = R_I*sin(interval)-0.25;
    plot(x1,x2);
    hold off
    title('L-shape domain');
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
    
    % Trefethen utilities
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