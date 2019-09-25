function  [k_m,sigma_min]=MFS_ITE_kite_isotropic(m,init,epsilon,Rrel_S,n)
% ROUTINE: Computes approximate ITE for Cakoni-ellipse with modified MFS
% INPUT:  
%   m = number of boundary collocation points (or equivalently source points)
%   init = initial guess for approximate ITE location
%   epsilon = ellipse perturbation parameter    
%   Rrel_S = scaling factor for scaling-translated sources wrt kite
%   n = refractive index

    format long
    
% Interior points parameters
    % Number of interior points
    m_I=20;
    % Radius of circular interior points
    R_I=0.2; % floating radius would be "0.5*(0.75-epsilon)"

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
    boundary_x=0.75*cos(phi)+epsilon*cos(2*phi);
    boundary_y=sin(phi);
    boundary=[boundary_x ; boundary_y];
    
    % Normal vectors along boundary
    normal_x=cos(phi)./sqrt(cos(phi).^2+(0.75*sin(phi)+2*epsilon*sin(2*phi)).^2);
    normal_y=(0.75*sin(phi)+2*epsilon*sin(2*phi))./sqrt(cos(phi).^2+(0.75*sin(phi)+2*epsilon*sin(2*phi)).^2);
    normals=[normal_x ; normal_y];
    
    % Source points 
    sources=[Rrel_S*boundary_x ; Rrel_S*boundary_y];
        
% Modified MFS main routine
    f_tmp=@(kappa) (f(kappa,m,m_I,interior,boundary,sources,normals,n)); 
    [k_m,sigma_min]=fminsearch(f_tmp,init,optimset('TolX',1e-16));
    %[k_m,sigma_min]=fminbnd(f_tmp,init-0.01,init+0.01,optimset('TolX',1e-16));
    
%     figure(m)
%     fplot(f_tmp, [init-0.001 init+0.001]);

%   Plot domain with sampled scattering points  
    plot_flag=false;
    if (plot_flag==true)
        interval = linspace(0,2*pi);
        y1 = 0.75*cos(interval)+epsilon*cos(2*interval);
        y2 = sin(interval);
        figure(1)
        plot(y1,y2);
        axis equal;
        hold on
        n_sparse=m/6;
        d_sparse=2*pi/(4*n_sparse);
        d_dense=d_sparse/2.;
        gap=d_dense/2;
        distr=[0:d_dense:(pi/2-gap) (pi/2):d_sparse:(pi-gap) pi:d_dense:(1.5*pi-gap) (1.5*pi):d_sparse:(2*pi-gap)];
        phi=distr+ones*pi/4;
        x1 = R_I*cos(interval)+xshift;
        x2 = R_I*sin(interval);
        z1 = R_S*cos(interval);
        z2 = R_S*sin(interval);
        plot(x1,x2);
        hold on
        plot(z1,z2);
        hold on
        quiver(0.75*cos(phi)+epsilon*cos(2*phi),sin(phi),cos(phi)./sqrt(cos(phi).^2+(0.75*sin(phi)+2*epsilon*sin(2*phi)).^2)...
                                               ,(0.75*sin(phi)+2*epsilon*sin(2*phi))./sqrt(cos(phi).^2+(0.75*sin(phi)+2*epsilon*sin(2*phi)).^2));
        title('Cakoni ellipse & normal vectors');                                   
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
    
    % Setup up the interior points of interior domain 
    M_I_i=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        M_I_i(i,j)=kernelHelmholtz(pts_I(:,i),pts_S(:,j),lambi);
      end
    end
    
    % Setup up the interior points of exterior domain 
    M_I_e=zeros(m_I,m);
    for i=1:m_I
      for j=1:m
        M_I_e(i,j)=kernelHelmholtz(pts_I(:,i),pts_S(:,j),sqrt(n)*lambi);
      end
    end
    
    % Trefethen
    A=[M_D_i  M_D_e ; M_N_i  M_N_e ; M_I_i  zeros(m_I,m) ; zeros(m_I,m) M_I_e];
    [Q,~]=qr(A,0);
    s=min(svd(Q(1:(2*m),:)));
end

function z=kernelHelmholtz(x,y,k)
    r=norm(x-y,2);
    z=besselh(0,1,k*r);
end

function z=nabla_kernelHelmholtz(x,y,k) 
    r=norm(x-y,2);
    z=-k*besselh(1,1,k*r)*(x-y)'/r;
end
