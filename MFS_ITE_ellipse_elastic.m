function  [w_m,sigma_min]=MFS_ITE_ellipse_elastic(m,init,semi_x,semi_y,Rrel_S,mu,lambda,rho)

% ROUTINE: Computes approximate elastic ITEs for ellipse scatterer via modified MFS
% INPUT:  
%   m = number of boundary collocation points (or equivalently source points)
%   init = initial guess for approximate ITE location
%   semi_x = semi axis of ellipse in x direction
%   semi_y = semi axis of ellipse in y direction    
%   Rrel_S = scaling factor of radius for circular source point radius wrt semi major axis
%   mu = Lamé parameter mu
%   lambda = Lamé parameter lambda
%   rho = density of scattering medium

    format long
    
% Interior & source point parameters
    % Number of interior points
    m_I=10; 
    % Scaling factor wrt minor axis for radius R_I of circular interior points
    Rrel_I=0.5;
    % Absolute radii for inner and source points
    R_I=Rrel_I*min(semi_x,semi_y);
    R_S=Rrel_S*max(semi_x,semi_y);
    
% Generation of computational points arrays 
    % Interior points 
    phi_I=linspace(0,2*pi,m_I+1); 
    phi_I=phi_I(1:end-1);
    interior_x=cos(phi_I);
    interior_y=sin(phi_I);
    interior=[R_I*interior_x ; R_I*interior_y];  
    
    % Collocation points
    phi=linspace(0,2*pi,m+1);
    phi=phi(1:end-1); 
    collocation=[semi_x*cos(phi); semi_y*sin(phi)];
    % Normal vectors along collocation boundary
    normals_x=semi_y*cos(phi)./sqrt(semi_y^2*cos(phi).^2+semi_x^2*sin(phi).^2);
    normals_y=semi_x*sin(phi)./sqrt(semi_y^2*cos(phi).^2+semi_x^2*sin(phi).^2);
    normals=[normals_x ; normals_y];
    % Source points (scaled version of collocation boundary)
    sources=[R_S*cos(phi) ; R_S*sin(phi)];

% Modified MFS main routine
    f_tmp=@(omega) (f(omega,m,m_I,interior,collocation,sources,normals,mu,lambda,rho)); 
    [w_m,sigma_min]=fminsearch(f_tmp,init,optimset('TolX',1e-16));%,'TolFun',1e-16));
    %[w_m,sigma_min]=fminbnd(f_tmp,init,init+0.5,optimset('TolX',1e-16));  
          
%%   Plot scattering domain, boundary normals, interior and source points
     %figure(m)
     %fplot(f_tmp, [0 2]);

end


function s=f(w,m,m_I,pts_I,pts_C,pts_S,normals,mu,lambda,rho)

    % Exception since besselh(1,1,w) is singular at w=0
    if(w==0)
        s=0;
        return
    end
 
    % Setup up Dirichlet boundary for u
    D_u=zeros(2*m,2*m);
    for ii=1:m
      for jj=1:m
        D_u(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=fundamental(pts_C(:,ii),pts_S(:,jj),w,lambda,mu);
      end
    end      
    
    % Setup up Dirichlet boundary for v
    D_v=zeros(2*m,2*m);
    for ii=1:m
      for jj=1:m
        D_v(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=fundamental(pts_C(:,ii),pts_S(:,jj),w*sqrt(rho),lambda,mu);
      end
    end  
    
    % Setup up Neumann boundary for u
    N_u=zeros(2*m,2*m);
    for ii=1:m
      for jj=1:m
        N_u(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=stresstensor(pts_C(:,ii),pts_S(:,jj),w,lambda,mu,normals(:,ii));
      end
    end      
    
    % Setup up Neumann boundary for v
    N_v=zeros(2*m,2*m);
    for ii=1:m
      for jj=1:m
        N_v(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=stresstensor(pts_C(:,ii),pts_S(:,jj),sqrt(rho)*w,lambda,mu,normals(:,ii));
      end
    end      
    
    % Setup up interior samples of u
    I_u=zeros(2*m_I,2*m);
    for ii=1:m_I
      for jj=1:m
        I_u(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=fundamental(pts_I(:,ii),pts_S(:,jj),w,lambda,mu);
      end
    end
    
    % Setup up interior samples of v
    I_v=zeros(2*m_I,2*m);
    for ii=1:m_I
      for jj=1:m
        I_v(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=fundamental(pts_I(:,ii),pts_S(:,jj),w*sqrt(rho),lambda,mu);
      end
    end
    
    % Trefethen QR utilities
    A=[D_u  D_v ; N_u  N_v ; I_u  zeros(2*m_I,2*m);  zeros(2*m_I,2*m) I_v];
    [Q,~]=qr(A,0);
    tmp=svd(Q(1:(4*m),:));
    s=tmp(4*m);  
end


function z=fundamental(x,y,w,lambda,mu)    

    comp1=x(1)-y(1);
    comp2=x(2)-y(2);
    r = sqrt(comp1^2+comp2^2);
    
    kp=w/sqrt(lambda+2*mu);
    ks=w/sqrt(mu);
    
    s1 = 1i/(4*mu)*besselh(0,1,ks*r)...
           -1i/(4*w^2*r)*(ks*besselh(1,1,ks*r)-kp*besselh(1,1,kp*r));
    s2 = 1i/(4*w^2)*( 2./r*(ks*besselh(1,1,ks*r)-kp*besselh(1,1,kp*r))...
           -(ks^2*besselh(0,1,ks*r)-kp^2*besselh(0,1,kp*r)));
    J=[comp1^2 comp1*comp2 ; comp1*comp2 comp2^2]/r^2;
    z=s1.*eye(2)+s2.*J;
end


function z=stresstensor(x,y,omega,lambda,mu,normals)

    a=x(1)-y(1);
    b=x(2)-y(2);
    r=sqrt(a*a+b*b);
    
    ks=omega/sqrt(mu);
    kp=omega/sqrt(lambda+2*mu);
    
    H1ks=besselh(1,1,ks*r);
    H0ks=besselh(0,1,ks*r);
    H1kp=besselh(1,1,kp*r);
    H0kp=besselh(0,1,kp*r);
    phi2=1i/(4*omega^2)*(2*ks./r.*H1ks-ks^2*H0ks...
                                -2*kp./r.*H1kp+kp^2*H0kp)...
                                -(lambda+mu)/(4*pi*mu*(lambda+2*mu));
    phi2p=1i/(4*omega^2)*(ks^3*H1ks-4*ks*H1ks./r.^2+2*ks^2*H0ks./r...
                                 -kp^3*H1kp+4*kp*H1kp./r.^2-2*kp^2*H0kp./r);
    phi1p=-1i/(4*mu)*H1ks*ks-1i/(4*omega^2)*(H0ks*ks^2./r-2*ks*H1ks./r.^2-H0kp*kp^2./r+2*kp*H1kp./r.^2)...
                                                      +(lambda+3*mu)./(4*pi*mu*r*(lambda+2*mu));
    J=[a^2 a*b; a*b b^2]/r^2;
    vec=[a;b];
    EYE=eye(2,2);
    U1=       lambda*normals*vec'+mu*vec*normals'+mu*normals'*vec*EYE;
    U2=(lambda+2*mu)*normals*vec'+mu*vec*normals'+mu*normals'*vec*(EYE-4*J);
    F=phi1p/r*U1+phi2p/r*U1*J+phi2/r^2*U2;
    
    T=((eye(2,2)+2*J*(lambda+mu)/mu)*(-(a*normals(1)+b*normals(2))/r^2)+(-[a;b]*normals'-(-[a;b]*normals')')/r^2)*mu/(2*pi*(lambda+2*mu));
    z=F+T;
end
