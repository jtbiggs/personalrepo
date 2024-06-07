%Dynamic Beam with clamped end condition
clear
%Initialization of Variables

a=0; %interval begin
b=1; %interval end

Nbasis=13;
N=Nbasis+2; %number of basis splines

ip = @(x) 0; %initial position t=0
iv = @(x) 50*sin(pi^2*x); %initial velocity
Cl=@(x) 0;
t0=0; %initial time
tf=5; %final time

gamma1=0.025; %Viscous Damping
gamma2=100; %Kelvin-Voight damping coefficients
I=1.734*10^(-7); %beam moment of interia
E=2*10^6; %young's modulus
rho=980; %density of beam
a=0.032; %cross sectional area
g=-9.81; %gravity
mb=1.927; %mass beam
l=b-a; %beam length

M=zeros(N); %mass matrix
K=zeros(N); %stiffness matrix
D=zeros(N); %damping matrix
F=zeros(N,1); %F vector - integral of f(0,x)*B_i over [a, b]
G=zeros(N,1);
Fbar=zeros(N,1);
Gbar=zeros(N,1);
IV=zeros(N,1);
IP=zeros(N,1);
%Generation of Basis functions
[Basis,D_Basis,DD_Basis,~,~,~,~,~,~,~,~,~,xarr,~]=cubicbspline2(a,b,Nbasis);

%Generation of mass and stiffness matrices
for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,rho*a*Basis(i,:).*Basis(j,:));
        K(i,j)=trapz(xarr,E*I*DD_Basis(i,:).*DD_Basis(j,:));
        D(i,j)=trapz(xarr,gamma1*Basis(i,:).*Basis(j,:))+...
            trapz(xarr,gamma2*I*DD_Basis(i,:).*DD_Basis(j,:));
    end
end

%Creation of matrix A in xdot(t)=Ax(t)
A=zeros(2*N);
A=[zeros(N) eye(N); -M\K -M\D];

%Related to initial condition - f(0,x) and g(0,x)

alpha0=zeros(N,1); %initial position

for j=1:N
    IP(j,:)=trapz(xarr,ip(xarr).*Basis(j,:)); %initial position
    IV(j,:)=trapz(xarr,iv(xarr).*Basis(j,:)); %initial velocity
    Gbar(j,:)=trapz(xarr,mb*g/l.*Basis(j,:));
    Fbar(j,:)=trapz(xarr,Cl(xarr).*Basis(j,:));
end

G=[zeros(N,1);M\Gbar];
F=[zeros(N,1);M\Fbar];
alpha0=inv(M)*IP;
beta0=inv(M)*IV; %initial velocity (g = 0)

y_init=zeros(2*N,1);
y_init=[alpha0; beta0]; %initial condition vector

tspan=t0:0.05:tf;
[t,y_yprime]=ode15s(@dbeam,tspan,y_init,[],A,G,F); %solves system of ordinary differential equations

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
        Sol(i,:)=Sol(i,:)+y(i,j).*Basis(j,:);
    end
end

mesh(xarr,t,Sol) %Use in presentation
axis([0 1 t0 tf -10 10])
title('Clamped Dynamic Beam With LQR Controller')
xlabel('length - s')
ylabel('time - t')
zlabel('Deflection - w(s,t) (m)')


%Exact Solution
%exact_sol=DynamicbeamExact(a,b,N,t);

%plot of numerical solution vs exact solution at t=1
% plot(xarr,Sol(length(t),:),'--o','MarkerIndices',1:500:length(xarr))
% hold on
% plot(xarr,exact_sol(length(t),:))
% title('Clamped Dynamic Beam Solutions at t=1')
% xlabel('x')
% ylabel('Beam Deflection - y(x)')
% legend({'Numerical Solution','Exact Solution'},'Location','northeast')
% axis([0 1 -1 1])

