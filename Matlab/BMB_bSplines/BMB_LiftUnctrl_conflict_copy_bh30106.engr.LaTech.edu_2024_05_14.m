
clear, clc

%Static Variables
%a=0; %bmb start
%b=2; %bmb end
t0=0;
tf=10;



wdot0 = @(x) -2; %initial vertical velocity at t=0


%beam properties
l = .6096;  %0.6096 m = 2 ft, length of both beams put together 
l1 = l/2; %m, 1 ft, splitting total length between two wings
l2 = l1;  %m, 1 ft, splitting total length between two wings
lm = .0508; %=2in  %.1016 m = 4 in, length of rigid mass
%so entire spatial domain is [0, l1+lm+l2]
total_l = l1+lm+l2;
%lendrbeam = l1+lm;




num_b_splines=9; 
N = num_b_splines; %used to generate N total splines over [a,b] with (N-1)/2 splines on each beam and one shared middle spline

[~,~,~,~,~,~,~,~,~,BMBBasis,D_BMBBasis,DD_BMBBasis,xarr,x_spacing] = cubicbspline2(0,l,num_b_splines); %spline generator fcn
%plot(xarr,BMBBasis)
%Stiffness and Mass Matrices
K=zeros(N);  %K(i,j)=int_a^b B''(i)*B''(j) dx
M=zeros(N);  %M(i,j)=int_a^b B(i)*B(j) dx
D=zeros(N);
w0=@(x) 0;
%w0= @(x) (x-total_l/2).^2; %initial position at t=0

%Physical Characteristics
%assuming 1.9247 kg total mass of whole beam
m = 1.927;            %kg, mass of point load
mass = 1.927;          %kg, distributed load
beamdepth = .127;        %.127 m = 5 in
beamheight = .0254;      %.0254 m = 1 in
a = beamdepth*beamheight;   %m^2, cross-sectional area of beam
E = 2e6; %1.03e10;          %N/m^2, Young's modulus
rho = 980;%198.1238; %1600;    %density of beam material
I = (beamdepth*beamheight^3)/12; %m^4, area moment of inertia
gamma1 = 0.025;%0.01;%0;       %viscous damping coefficient
gamma2 = 1e2;%1e4;%0;%0.001;   %Kelvin-Voight damping coefficient
g = 9.81;                      %m/s^2, gravity at sea level 

rhoA = rho*a;
EI = E*I;
gamma2I = gamma2*I;
dist_mg = mass*g;
mgbyl = dist_mg/l;

%rigid mass properties
massdepth = beamdepth;  %0.0762;      %m = 3 in
massheight = beamheight;  %0.0508;    %m = 2 in 
%kg m^2, mass moment of inertia for rigid mass
Iz = (1/12)*m*(massdepth^2+massheight^2);


%Lift Coefficients
u=25; %m/sec (velocity of vehicle)
c=0.127; %inches
rho_a = 1.2; %density of air at sea level
k1=0.2; %0
k2=1.6;
k3=2.1;
k4=0; %7.2
k5=1.47068;       %wind velocity
liftc = -0.5*rho_a*u^2*c/(l); %coefficients for lift function

W0=zeros(N,1);
Wdot0=zeros(N,1);
Bbar=zeros(N,1);
G1=zeros(N,1);
for i=1:N
    for j=1:N
        %mass matrix
        M(i,j)=trapz(xarr,rhoA*BMBBasis(i,:).*BMBBasis(j,:));+...
            m*BMBBasis(i,(length(xarr)-1)/2+1)*BMBBasis(j,(length(xarr)-1)/2+1)+...
            Iz*D_BMBBasis(i,(length(xarr)-1)/2+1)*D_BMBBasis(j,(length(xarr)-1)/2+1)';
        %stiffness matrix
        K(i,j)=trapz(xarr,EI*DD_BMBBasis(i,:).*DD_BMBBasis(j,:));
        %kelvin voight damping matrix
        D(i,j)=trapz(xarr,gamma1*BMBBasis(i,:).*BMBBasis(j,:))+...
            trapz(xarr,gamma2*I*DD_BMBBasis(i,:).*DD_BMBBasis(j,:));
    end
    %Control vector
    Bbar(i)=trapz(xarr,BMBBasis(i,:));
    %gravity vector
    G1(i)=trapz(xarr,mgbyl*BMBBasis(i,:));
end
%M((N+1)/2,(N+1)/2)=M((N+1)/2,(N+1)/2)+2*m*36;
%Weak form initial displacement and velocity
for i=1:N
    W0(i,:)=trapz(xarr,w0(xarr).*BMBBasis(i,:));
    Wdot0(i,:)=trapz(xarr,wdot0(xarr).*BMBBasis(i,:));
end

%Vectorize for use in ode solver
%A matrix - xdot=A*x(t)
A=zeros(2*N);
A=[zeros(N) eye(N); (-M)\K (-M)\D];
%G Matrix - [0;M^(-1)*G1]
G=zeros(2*N,1);
G(N+1:2*N,1)=(M\G1);

alpha0=zeros(N,1); %initial position vector
alpha0=M\W0;
beta0=zeros(N,1); 
beta0=M\Wdot0; 

%initial vector
%beta0(1:2:end)=0;

%LQR Control Dynamics
% Q=2*eye(2*N);
% R=2;
% B=zeros(2*N,1);
% B(N+1:2*N)=M\Bbar;

%initial guess for ARE
% X=eye(2*N); 
%solutions for steady state ARE problem
% P=newt(A,B,Q,inv(R),X); 
% K_lqr=inv(R)*B'*P; 
%K gain matrix and Ricatti equation solutions
%[K_lqr,P]=lqr(A,B,Q,R); 


y_init=zeros(2*N,1);
y_init=[alpha0; beta0]; %initial condition vector
tspan=t0:0.01:tf;

%test=liftfun(tspan,y_init,BMBBasis,xarr,A,G,k1,k2,k3,k4,k5,liftc,u,N,M);

[t,alphac]=ode15s(@liftfun,tspan,y_init,[],BMBBasis,xarr,A,G,k1,k2,k3,k4,k5,liftc,u,N,M); %solves system of ordinary differential equations

y=zeros(length(t),N); %y sol coef
yprime=zeros(length(t),N); %yprime sol coef
for i=1:N
y(:,i)=alphac(:,i);
end
for j=1:N
yprime(:,j)=alphac(:,N+j);
end

%Solution Matrix for [a,b] and [t0,tf]
Sol=zeros(length(t),length(xarr));
for i=1:length(t)
    for j=1:N
        Sol(i,:)=Sol(i,:)+y(i,j).*BMBBasis(j,:);
    end
end

%Plot of BMB system as a function of s and t

figure('Name','BMB System')
[xmass,Solmass,mp,xgap]=massplot(0,l,x_spacing,lm,Sol,N,t);

%mesh(xmass(1:mp),t,Solmass(:,1:mp))
%hold on
%grid on
%mesh(xmass(mp+1:mp+xgap-1),t,Solmass(:,(mp+1:mp+xgap-1)),'EdgeColor','#b4b4b4')
%mesh(xmass(mp+xgap:length(xmass)),t,Solmass(:,mp+xgap:length(xmass)))
% 
axis([0 total_l  t0 tf min(Sol,[],'all') max(Sol,[],'all')])
title('Uncontrolled BMB Model with Lift')
xlabel('length - x')
zlabel('Deflection - w(x,t)')
ylabel('time - t')

%Plot of control effort by LQR controller over time
% figure('Name','Control Effort')
% u_eff=-K_lqr*y_yprime';
% plot(t,u_eff)
% title('Control Effort vs. Time')
% ylabel('Control Effort - u(t)')
% xlabel('time - t')

%Plot of BMB system at t=5
% figure('Name','Clamped Dynamic BMB Model with LQR Controller t=5')
% plot(xmass,Solmass(101,:))
% axis([0 2 -1 1])
% title('Clamped Dynamic BMB Model with LQR Controller t=5')
% xlabel('Length (m) - s')
% ylabel('Deflection (m) - w(s,t)')

