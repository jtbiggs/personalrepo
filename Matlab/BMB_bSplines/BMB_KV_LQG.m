%Dynamic BMB model with clamped left side
clc

%Static Variables
a=0; %bmb start
b=2; %bmb end
t0=0;
tf=10;
mass_size=0.25;
f=@(x) -0.5*x; %initial condition at t=0
num_b_splines=25; 
N = num_b_splines-2; %used to generate N total splines over [a,b] with (N-1)/2 splines on each beam and one shared middle spline

[~,~,~,~,~,~,ClampedBMBBasis,D_ClampedBMBBasis,DD_ClampedBMBBasis,xarr,x_spacing] = cubicbspline(a,b,num_b_splines); %spline generator fcn

%Stiffness and Mass Matrices
K=zeros(N);  %K(i,j)=int_a^b B''(i)*B''(j) dx
M=zeros(N);  %M(i,j)=int_a^b B(i)*B(j) dx
D=zeros(N);

gamma1=0.025; %Kelvin-Voight damping coefficients
gamma2=100;
I=1.734*10^(-7);

F=zeros(N,1);
Bbar=zeros(N,1);
for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,ClampedBMBBasis(i,:).*ClampedBMBBasis(j,:));
        K(i,j)=trapz(xarr,DD_ClampedBMBBasis(i,:).*DD_ClampedBMBBasis(j,:));
        D(i,j)=trapz(xarr,gamma1*ClampedBMBBasis(i,:).*ClampedBMBBasis(j,:))+...
            trapz(xarr,gamma2*I*DD_ClampedBMBBasis(i,:).*DD_ClampedBMBBasis(j,:));
    end
    Bbar(i)=trapz(xarr,ClampedBMBBasis(i,:));
end

for i=1:N
    F(i,:)=trapz(xarr,f(xarr).*ClampedBMBBasis(i,:));
end

%A matrix - xdot=A*x(t)  System Dynamics
A=zeros(2*N);
A=[zeros(N) eye(N); (-M)\K (-M)\D];

alpha0=zeros(N,1); %initial position vector
alpha0=M\F;
beta0=zeros(N,1); %initial velocity (g = 0)

%LQG Controller____________________________________________________________

%LQG Control Dynamics
C=2*[eye(N) zeros(N); zeros(N) zeros(N)];%obervability matrix %setting C=I reverts to LQR problem
Q=10*(C*C');
R=2;
B=zeros(2*N,1);
B(N+1:2*N)=M\Bbar;
x0=zeros(N,1);
x0=[alpha0;beta0]; %initial condition vector

X=eye(2*N); %initial guess for ARE

pie_lqg=newt(A,B,Q,inv(R),X);  %steady state ARE solution
fbK_lqg=R\B'*pie_lqg;           %feedback K gain matrix 
pea_lqg = are(A',C'*C,B*B');
F_lqg = pea_lqg*C';
Ac_lqg = A-B*fbK_lqg-F_lqg*C;
A_lqg = [A  -B*fbK_lqg;F_lqg*C  Ac_lqg];

xspan=0:b/(4*N-1):b;
xend = zeros(1,2*N); %regulator problem if all zeros.. otherwise tracking problem
x_end = [xend xend]';

tspan = t0:0.1:tf;
%Get little b for feedforward term
backwards_time = tf:-0.05:0;
b_end = zeros(1,length(B'));
[tb, bb] = ode15s('findb',backwards_time,b_end,[],A,B,inv(R),pie_lqg,x_end);
bb=bb(end,:)'; %Steady state solution
u_fw = B*inv(R)*B'*bb; %feedforward term


d0=1;
x0est=[x0; 0.75*x0];
%%%____________________________________________________________________________

[t,y_yprime]=ode45(@(t,x) dbeam_lqg(t,x,A_lqg,B,fbK_lqg,u_fw,F,d0,x_end),tspan,x0est); %solves system of ordinary differential equations

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
grid on
mesh(xmass(mp+1:mp+xgap-1),t,Solmass(:,(mp+1:mp+xgap-1)),'EdgeColor','#b4b4b4')
mesh(xmass(mp+xgap:length(xmass)),t,Solmass(:,mp+xgap:length(xmass)))

axis([0 2+mass_size  t0 tf -3 3])
title('Clamped Dynamic BMB Model with LQG Controller')
xlabel('length - x')
zlabel('Deflection - w(x,t)')
ylabel('time - t')
