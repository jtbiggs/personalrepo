%Dynamic BMB model with clamped left side
clc

%Static Variables
a=0; %bmb start
b=2; %bmb end
t0=0;
tf=15;
mass_size=0.25;
f=@(x) -0.5*x; %initial condition at t=0
num_b_splines=19; 
N = num_b_splines-2; %used to generate N total splines over [a,b] with (N-1)/2 splines on each beam and one shared middle spline

[~,~,~,~,~,~,ClampedBMBBasis,D_ClampedBMBBasis,DD_ClampedBMBBasis,xarr,x_spacing] = cubicbspline(a,b,num_b_splines); %spline generator fcn

%Stiffness and Mass Matrices
K=zeros(N);  %K(i,j)=int_a^b B''(i)*B''(j) dx
M=zeros(N);  %M(i,j)=int_a^b B(i)*B(j) dx
F=zeros(N,1);
Bbar=zeros(N,1);
for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,ClampedBMBBasis(i,:).*ClampedBMBBasis(j,:));
        K(i,j)=trapz(xarr,DD_ClampedBMBBasis(i,:).*DD_ClampedBMBBasis(j,:));
    end
    Bbar(i)=trapz(xarr,ClampedBMBBasis(i,:));
end

for i=1:N
    F(i,:)=trapz(xarr,f(xarr).*ClampedBMBBasis(i,:));
end

%A matrix - xdot=A*x(t)
A=zeros(2*N);
A=[zeros(N) eye(N); (-M)\K zeros(N)];

alpha0=zeros(N,1); %initial position vector
alpha0=M\F;
beta0=zeros(N,1); %initial velocity (g = 0)

%LQR Control Dynamics
Q=eye(2*N);
R=1;
B=zeros(2*N,1);
B(N+1:2*N)=M\Bbar;
[K_lqr,P]=lqr(A,B,Q,R); %Generates K gain matrix and Ricatti equation solutions


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
        Sol(i,:)=Sol(i,:)+y(i,j).*ClampedBMBBasis(j,:);
    end
end

%Plot of BMB system as a function of s and t
figure('Name','BMB System')
[xmass,Solmass,mp,xgap]=massplot(a,b,x_spacing,mass_size,Sol);

mesh(xmass(1:mp),t,Solmass(:,1:mp))
hold on
mesh(xmass(mp+1:mp+xgap-1),t,Solmass(:,(mp+1:mp+xgap-1)))
mesh(xmass(mp+xgap:length(x_mass)),t,Solmass(:,mp+xgap:length(x_mass)))


axis([0 2  0 15 -3 3])
title('Clamped Dynamic BMB Model with LQR Controller')
xlabel('length (m) - s')
zlabel('Deflection (m) - w(s,t)')
ylabel('time (s) - t')
grid

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
