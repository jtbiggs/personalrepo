%Dynamic Beam with clamped end condition
clc
%Initialization of Variables

a=0; %interval begin
b=1; %interval end
N=23; %number of basis splines

f = @(x) -0.5*x; %initial position t=0
g = 0; %initial velocity
t0=0; %initial time
tf=5; %final time

gamma1=0.025; %Kelvin-Voight damping coefficients
gamma2=100;
I=1.734*10^(-7);

M=zeros(N); %mass matrix
K=zeros(N); %stiffness matrix
D=zeros(N); %damping matrix
F=zeros(N,1); %F vector - integral of f(0,x)*B_i over [a, b]

%Generation of Basis functions
[Basis,D_Basis,DD_Basis,ClampedBasis,D_ClampedBasis,DD_ClampedBasis,~,~,~,xarr]=cubicbspline(a,b,N);

%Generation of mass and stiffness matrices
for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,ClampedBasis(i,:).*ClampedBasis(j,:));
        K(i,j)=trapz(xarr,DD_ClampedBasis(i,:).*DD_ClampedBasis(j,:));
        D(i,j)=trapz(xarr,gamma1*ClampedBasis(i,:).*ClampedBasis(j,:))+...
            trapz(xarr,gamma2*I*DD_ClampedBasis(i,:).*DD_ClampedBasis(j,:))+...
            trapz(xarr,xarr.*ClampedBasis(i,:).*ClampedBasis(j,:));

    end
end

%Creation of matrix A in xdot(t)=Ax(t)
A=zeros(2*N);
A=[zeros(N) eye(N); -inv(M)*K -inv(M)*D];

%Related to initial condition - f(0,x) and g(0,x)

alpha0=zeros(N,1); %initial position
for i=1:N
    F(i,:)=trapz(xarr,f(xarr).*ClampedBasis(i,:));
end
alpha0=inv(M)*F;

beta0=zeros(N,1); %initial velocity (g = 0)

y_init=zeros(2*N,1);
y_init=[alpha0; beta0]; %initial condition vector

tspan=t0:0.01:tf;
[t,y_yprime]=ode15s(@dbeam,tspan,y_init,[],A); %solves system of ordinary differential equations

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
        Sol(i,:)=Sol(i,:)+y(i,j).*ClampedBasis(j,:);
    end
end

mesh(xarr,t,Sol) %Use in presentation
axis([0 1 t0 tf -1 1])
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

