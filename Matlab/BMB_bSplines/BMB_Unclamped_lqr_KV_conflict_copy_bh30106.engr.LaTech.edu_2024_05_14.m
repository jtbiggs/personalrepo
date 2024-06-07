%Dynamic BMB model with clamped left side
clear

%Static Variables
a=0; %bmb start
b=2; %bmb end
t0=0;
tf=5; %simulation end time


f=@(x) (x-1).^2;%initial displacement at t=0
g= @(x) 0; %initial velocity
num_b_splines=21; %number of nodes
N = num_b_splines; %used to generate N total splines over [a,b] with (N-1)/2 splines on each beam and one shared middle spline
sys_dof=2*N;

[~,~,~,~,~,~,~,~,~,BMBBasis,D_BMBBasis,DD_BMBBasis,xarr,x_spacing] = cubicbspline2(a,b,num_b_splines); %spline generator fcn

%Stiffness and Mass Matrices
K=zeros(N);  %K(i,j)=int_a^b B''(i)*B''(j) dx
M=zeros(N);  %M(i,j)=int_a^b B(i)*B(j) dx
init_pos=zeros(N,1);
init_vel=zeros(N,1);
D=zeros(N);
G=zeros(N,1);
Bbar=zeros(N,1);

%%%Beam Properties
mass_size=.0508;
l1=(b-a)/2; %length of left beam
l2=l1; %length of right beam
l_beams=l1+l2; %length of beams together
total_l=l1+l2+mass_size; %total length of system

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
num_of_gauss = 6;

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

%lift function parameters
u=25; %m/sec (velocity of vehicle)
c=0.127; %inches
rho_a = 1.2; %density of air at sea level
k1=0.2; %0
k2=1.6;
k3=2.1;
k4=0; %7.2
k5=1.47068;       %wind velocity
liftc = -0.5*rho_a*u^2*c/(l_beams); %coefficients for lift function



for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,rhoA*BMBBasis(i,:).*BMBBasis(j,:));
        K(i,j)=trapz(xarr,EI*DD_BMBBasis(i,:).*DD_BMBBasis(j,:));
        D(i,j)=trapz(xarr,gamma1*BMBBasis(i,:).*BMBBasis(j,:))+...
            trapz(xarr,gamma2*I*DD_BMBBasis(i,:).*DD_BMBBasis(j,:));
        
    end
    Bbar(i)=trapz(xarr,BMBBasis(i,:));
    G(i)=trapz(xarr,(m*g)/l_total*BMBBasis(i,:));
end

for i=1:N
    init_pos(i,:)=trapz(xarr,f(xarr).*BMBBasis(i,:));
    init_vel(i,:)=trapz(xarr,g(xarr).*BMBBasis(i,:));
end

%A matrix - xdot=A*x(t)
A=zeros(sys_dof);
A=[zeros(N) eye(N); (-M)\K (-M)\D];

alpha0=zeros(N,1); %initial position vector
alpha0=M\init_pos;

beta0=zeros(N,1); %initial velocity (g = 0)
beta0=M\init_vel;
%LQR Control Dynamics
Q=5*eye(sys_dof);
R=1;
B=zeros(sys_dof,1);
B(N+1:sys_dof)=M\Bbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1st Integral for Lift Function
%The rest will be assembled by 

F1=liftc*k1*trapz(xarr,BMBBasis); %first integral for lift function








X=eye(sys_dof); %initial guess for ARE
P=newt(A,B,Q,inv(R),X); %solutions for steady state ARE problem
K_lqr=inv(R)*B'*P; 

y_init=zeros(2*N,1);
y_init=[alpha0; beta0]; %initial condition vector
tspan=t0:0.05:tf;
[t,y_yprime]=ode15s(@dbeam_lqr,tspan,y_init,[],A,K_lqr,B); %solves system of ordinary differential equations

y=zeros(length(t),N); %y sol coef
yprime=zeros(length(t),N); %yprime sol coef
for i=1:N
y(:,i)=y_yprime(:,i);
end
for j=1:N
yprime(:,j)=y_yprime(:,N+j);
end

%Solution Matrix for [a,b] and [t0,tf]
Sol=zeros(length(t),length(xarr));
for i=1:length(t)
    for j=1:N
        Sol(i,:)=Sol(i,:)+y(i,j).*BMBBasis(j,:);
    end
end

figure('Name','BMB Initial Position')
% plot(xarr,Sol(1,:))
% hold on
% plot(xarr,f(xarr));


[xmass,Solmass,mp,xgap]=massplot(a,b,x_spacing,mass_size,Sol);

mesh(xmass(1:mp),t,Solmass(:,1:mp))
hold on
grid on
mesh(xmass(mp+1:mp+xgap-1),t,Solmass(:,(mp+1:mp+xgap-1)),'EdgeColor','#b4b4b4')
mesh(xmass(mp+xgap:length(xmass)),t,Solmass(:,mp+xgap:length(xmass)))

axis([a b+mass_size  t0 tf -1 1])
title('Uncontrolled Clamped Dynamic BMB Model')
xlabel('length - s')
zlabel('Deflection - w(s,t)')
ylabel('time - t')


