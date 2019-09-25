function  [k_m,sigma_min]=MFS_ITE_polygon_isotropic(m,init,p,Rrel_S,n)
% ROUTINE: Computes approximate ITE for regular p_Polygon of unit length with modified MFS
% INPUT:  
%   m = number of boundary collocation points (or equivalently source points)
%   init = initial guess for approximate ITE location
%   p = number of vertices of polygon    
%   Rrel_S = scaling factor for scaling-translated sources wrt outer polygon radius
%   n = refractive index

    format long;
    
    if(mod(m,2)==0 && m~=p)
        disp('   Choose m odd to avoid sample points at vertices!');
        return
    end

% Interior points parameters
    % Number of interior points
    m_I=20;
    % Scaling factor wrt inner circumference of polygon for radius R_I of circular interior points
    Rrel_I=0.5;
    % Absolute radius for inner points
    R_I=Rrel_I*0.5/tan(2*pi/(2*p));
    
% Generation of computational points 
    % Interior points 
    phi_I=linspace(0,2*pi,m_I+1);
    phi_I=phi_I(1:end-1);
    interior_x=cos(phi_I);
    interior_y=sin(phi_I);
    interior=[R_I*interior_x ; R_I*interior_y]; 
    
    % Boundary collocation points
        % midpoint to midedge
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
            boundary_x(1,ii+1)=R*cos(2*pi*ii/m);
            boundary_y(1,ii+1)=R*sin(2*pi*ii/m);
            normal_x(ii+1)=cos(section*2*pi/p);
            normal_y(ii+1)=sin(section*2*pi/p);
        end
    boundary=[boundary_x ; boundary_y];
    
    phi=linspace(0,2*pi,m+1);
    phi=phi(1:end-1);
    sources=[Rrel_S*MV*cos(phi) ; Rrel_S*MV*sin(phi)];
    %sources=[Rrel_S*boundary_x ; Rrel_S*boundary_y]; % source boundary scaled version of collocation boundary
    
    normals=[normal_x ; normal_y];
    
    % Modified MFS main routine
    f_tmp=@(kappa) (f(kappa,m,m_I,interior,boundary,sources,normals,n)); 
    [k_m,sigma_min]=fminsearch(f_tmp,init,optimset('TolX',1e-16));
    %[k_m,sigma_min]=fminbnd(f_tmp,init-0.01,init+0.01,optimset('TolX',1e-16));

end

% Computes minimal singular value for arbitrary wavenumber lambi
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
    %cond(Q(1:(2*m),:))
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