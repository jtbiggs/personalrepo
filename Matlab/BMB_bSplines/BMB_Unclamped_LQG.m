%Dynamic BMB model with clamped left side
clear

%Static Variables
a=0; %bmb start
b=2; %bmb end
t0=0;
tf=5;
mass_size=0.25;

f=@(x) (x-1).^2;%initial condition at t=0
g= @(x) 0;
num_b_splines=21; 
N = num_b_splines; %used to generate N total splines over [a,b] with (N-1)/2 splines on each beam and one shared middle spline

[~,~,~,~,~,~,~,~,~,BMBBasis,D_BMBBasis,DD_BMBBasis,xarr,x_spacing] = cubicbspline2(a,b,num_b_splines); %spline generator fcn

%Stiffness and Mass Matrices
K=zeros(N);  %K(i,j)=int_a^b B''(i)*B''(j) dx
M=zeros(N);  %M(i,j)=int_a^b B(i)*B(j) dx
F=zeros(N,1);
G=zeros(N,1);
D=zeros(N);
Bbar=zeros(N,1);
gamma1=0.025; %Kelvin-Voight damping coefficients
gamma2=100;
I=1.734*10^(-7);

for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,BMBBasis(i,:).*BMBBasis(j,:));
        K(i,j)=trapz(xarr,DD_BMBBasis(i,:).*DD_BMBBasis(j,:));
        D(i,j)=trapz(xarr,gamma1*BMBBasis(i,:).*BMBBasis(j,:))+...
            trapz(xarr,gamma2*I*DD_BMBBasis(i,:).*DD_BMBBasis(j,:));
    end
    Bbar(i)=trapz(xarr,BMBBasis(i,:));
end

for i=1:N
    F(i,:)=trapz(xarr,f(xarr).*BMBBasis(i,:));
    G(i,:)=trapz(xarr,g(xarr).*BMBBasis(i,:));
end

%A matrix - xdot=A*x(t)
A=zeros(2*N);
A=[zeros(N) eye(N); (-M)\K (-M)\D];

alpha0=zeros(N,1); %initial position vector
alpha0=M\F;

beta0=zeros(N,1); %initial velocity (g = 0)
beta0=M\G;
%LQG Controller____________________________________________________________

%LQG Control Dynamics
C=10*[eye(N) zeros(N); zeros(N) zeros(N)];%obervability matrix %setting C=I reverts to LQR problem
Q=C*C';
R=1;
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

tspan = t0:0.025:tf;
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
title('Unclamped BMB Model LQG')
xlabel('length - s')
zlabel('Deflection - w(t,s)')
ylabel('time - t')

