%Dynamic Beam with clamped end condition
clc
%Initialization of Variables

a=0; %interval begin
b=1; %interval end
N=21; %number of basis splines

f = @(x) -0.5*x; %initial position t=0
g = 0; %initial velocity
t0=0; %initial time
tf=5; %final time

gamma1=0.025; %Kelvin-Voigt damping coefficients
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
        M(i,j)=simps(xarr,ClampedBasis(i,:).*ClampedBasis(j,:));
        K(i,j)=simps(xarr,DD_ClampedBasis(i,:).*DD_ClampedBasis(j,:));
        %D(i,j)=trapz(xarr,gamma1*ClampedBasis(i,:).*ClampedBasis(j,:))+...
            %trapz(xarr,gamma2*I*DD_ClampedBasis(i,:).*DD_ClampedBasis(j,:))+...
            %trapz(xarr,xarr.*ClampedBasis(i,:).*ClampedBasis(j,:));

    end
end

%Creation of matrix A in xdot(t)=Ax(t)
A=zeros(2*N);
A=[zeros(N) eye(N); -M\K zeros(N)];

%Related to initial condition - f(0,x) and g(0,x)

alpha0=zeros(N,1); %initial position
for i=1:N
    F(i,:)=simps(xarr,f(xarr).*ClampedBasis(i,:));
end
alpha0=M\F;

beta0=zeros(N,1); %initial velocity (g = 0)

y_init=zeros(2*N,1);
y_init=[alpha0; beta0]; %initial condition vector

tspan=t0:0.05:tf;
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
for i=1:length(t)
    for j=1:N
        Sol(i,:)=Sol(i,:)+y(i,j).*ClampedBasis(j,:);
    end
end
%Exact Solution
exact_sol=DynamicbeamExact(a,b,N,t);

%plot of numerical solution vs exact solution at t=5
plot(xarr,Sol(101,:),'--o','MarkerIndices',1:ceil(length(xarr)/N):length(xarr))
hold on
plot(xarr,exact_sol(101,:))
title('Uncontrolled PDE Beam Solutions at t=5')
xlabel('Length - x')
ylabel('Beam Deflection - w(x)')
legend({'Numerical Solution','Exact Solution'},'Location','northeast')
axis([a b -1 1])

% figure(1)
% mesh(xarr,t,Sol);
% title('Uncontrolled PDE Beam Numerical Solution')
% xlabel('Length - x')
% zlabel('Deflection - w(t,x)')
% ylabel('Time - t')
%axis([0 1 t0 tf -1 1])
% 
% figure(2)
% mesh(xarr,t,exact_sol);
% title('Uncontrolled PDE Beam Exact Solution')
% xlabel('Length - x')
% zlabel('Deflection - w(t,x)')
% ylabel('Time - t')
% axis([0 1 t0 tf -1 1])


