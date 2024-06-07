%Dynamic BMB model with clamped left side
clear

%Static Variables
a=0; %bmb start
b=2; %bmb end
t0=0;
tf=5;
mass_size=0.25;
m=2; %mass in kg

f=@(x) 0.5*(x-1).^2; %initial condition at t=0
g=@(x) 0;
gamma1=0.025; %Kelvin-Voigt damping coefficients
gamma2=100;
I=1.734*10^(-7);

num_b_splines=21;
N=num_b_splines;

[~,~,~,~,~,~,~,~,~,BMBBasis,D_BMBBasis,DD_BMBBasis,xarr,x_spacing] = cubicbspline2(a,b,num_b_splines); %spline generator fcn

%Stiffness and Mass Matrices
K=zeros(N);  %K(i,j)=int_a^b B''(i)*B''(j) dx
M=zeros(N);  %M(i,j)=int_a^b B(i)*B(j) dx
D=zeros(N);
F=zeros(N,1);
G=zeros(N,1);
Bbar=zeros(N,1);
for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,BMBBasis(i,:).*BMBBasis(j,:));
        K(i,j)=trapz(xarr,DD_BMBBasis(i,:).*DD_BMBBasis(j,:));
        D(i,j)=trapz(xarr,gamma1*BMBBasis(i,:).*BMBBasis(j,:))+...
              trapz(xarr,gamma2*I*DD_BMBBasis(i,:).*DD_BMBBasis(j,:))+...
              trapz(xarr,xarr.*BMBBasis(i,:).*BMBBasis(j,:));
    end
    Bbar(i)=trapz(xarr,BMBBasis(i,:));
end
for i=1:N
    F(i,:)=trapz(xarr,f(xarr).*BMBBasis(i,:));
    G(i,:)=trapz(xarr,g(xarr).*D_BMBBasis(i,:));
end

%A matrix - xdot=A*x(t)
A=zeros(2*N);
A=[zeros(N) eye(N); (-M)\K zeros(N)];
alpha0=zeros(N,1); %initial position vector
alpha0=M\F;

beta0=zeros(N,1); %initial velocity (g = 0)
beta0=M\G;
y_init=zeros(2*N,1);
y_init=[alpha0; beta0]; %initial condition vector
tspan=t0:0.025:tf;
disp('running ode45')
[t,y_yprime]=ode45(@dbeam,tspan,y_init,[],A); %solves system of ordinary differential equations

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
D_Sol=zeros(length(t),length(xarr));
for i=1:length(t)
    for j=1:N
        Sol(i,:)=Sol(i,:)+y(i,j).*BMBBasis(j,:);
        D_Sol(i,:)=D_Sol(i,:)+yprime(i,j).*BMBBasis(j,:);
    end
end

figure('Name','BMB System')
[xmass,Solmass,mp,xgap]=massplot(a,b,x_spacing,mass_size,Sol);

mesh(xmass(1:mp),t,Solmass(:,1:mp))
hold on
grid on
mesh(xmass(mp+1:mp+xgap-1),t,Solmass(:,(mp+1:mp+xgap-1)),'EdgeColor','#b4b4b4')
mesh(xmass(mp+xgap:length(xmass)),t,Solmass(:,mp+xgap:length(xmass)))

axis([a b+mass_size  t0 tf -1 1])
title('Uncontrolled BMB Model')
xlabel('length - s')
zlabel('Deflection - w(t,s)')
ylabel('time - t')




