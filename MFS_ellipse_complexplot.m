function  MFS_ellipse_complexplot(m,semi_x,semi_y,Rrel_S,n,a_xx,a_yy,a_xy)
%Plot contour plot in some complex plane interval for anisotropic acoustic ITP
   % option = optional flag for computation by either standard(false)or modified MFS(true)
       option = true;  % false works without Beyn reformulation pretty bad!!!!!!!!!!
      
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
    N=100;
    x_max=3.1;
    shift=1.7;
    X=[shift:x_max/N:(shift+x_max)];
    y_max=1.5;
    Y=[-y_max:y_max/N:y_max];
    Z=zeros(2*N+1,N+1);
    f_tmp=@(omega) (f(omega,m,m_I,interior,boundary,sources,normals,n,a_xx,a_yy,a_xy,option)); 
     for ix=0:1:N
         for iy=-N:1:N 
         Z(iy+1+N,ix+1)= f_tmp(shift+ix*(x_max/N)+1i*iy*(y_max/N));
         end
     end
    figure('Name','ellipse','NumberTitle','off');
    contour(X,Y,Z,0:0.03:0.3) 
    %title(['{\it m}=',num2str(m)]);
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;
    axis equal;
    xlabel('Re(k)');
    ylabel('Im(k)');
    hold off         
end

function s=f(k,m,m_I,pts_I,pts_B,pts_S,normals,n,a_xx,a_yy,a_xy,option_flag)

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
    
    A=[D_v  D_w ; N_v  N_w ; I_v  zeros(m_I,m) ; I_vx  zeros(m_I,m); I_vy  zeros(m_I,m); zeros(m_I,m) I_w; zeros(m_I,m) I_wx; zeros(m_I,m) I_wy];
    if (option_flag==true)
    % Trefethen QR utilities for modified MFS
    [Q,~]=qr(A,0);
    tmp=svd(Q(1:(2*m),:));
    s=tmp(2*m); 
    else
    % standard MFS
    tmp=svd(A(1:(2*m),:));
    s=tmp(2*m);  
    end
end

function z=kernelHelmholtz(x,y,k)
    r=norm(x-y,2);
    z=besselh(0,1,k*r);
end

function z=nabla_kernelHelmholtz(x,y,k) 
    r=norm(x-y,2);
    z=-k*besselh(1,1,k*r)*(x-y)'/r;
end
