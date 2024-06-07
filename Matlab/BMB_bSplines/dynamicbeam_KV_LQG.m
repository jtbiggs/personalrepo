%Dynamic Beam with LQG controller
%Dynamic Beam with clamped end condition
clear all
%Initialization of Variables

a=0; %interval begin
b=1; %interval end
N=15; %number of basis splines
mass_size=0.25;

gamma1=0.025; %Kelvin-Voight damping coefficients
gamma2=100;
I=1.734*10^(-7);


f = @(x) -0.5*x; %initial position t=0
g = 0; %initial velocity
t0=0; %initial time
tf=7; %final time

M=zeros(N); %mass matrix
D=zeros(N); %damping matrix
K=zeros(N); %stiffness matrix
F=zeros(N,1); %F vector - integral of f(0,x)*B_i over [a, b]
Bbar=zeros(N,1);
b=1; %constant controller
%Generation of Basis functions
[Basis,D_Basis,DD_Basis,ClampedBasis,D_ClampedBasis,DD_ClampedBasis,~,~,~,xarr]=cubicbspline(a,b,N);

%Generation of mass and stiffness matrices
for i=1:N
    for j=1:N
        M(i,j)=trapz(xarr,ClampedBasis(i,:).*ClampedBasis(j,:));
        K(i,j)=trapz(xarr,DD_ClampedBasis(i,:).*DD_ClampedBasis(j,:));
        %D(i,j)=trapz(xarr,gamma1*ClampedBasis(i,:).*ClampedBasis(j,:))+...
            %trapz(xarr,gamma2*I*DD_ClampedBasis(i,:).*DD_ClampedBasis(j,:));

    end
    Bbar(i)=trapz(xarr,ClampedBasis(i,:));
end



%Related to initial condition - f(0,x) and g(0,x)

alpha0=zeros(2*N,1); %initial position
for i=1:N
    F(i,:)=trapz(xarr,f(xarr).*ClampedBasis(i,:));
end


alpha0=inv(M)*F;
beta0=zeros(N,1); %initial velocity (g = 0)
x0=zeros(N,1);
x0=[alpha0;beta0]; %initial condition vector
tspan=t0:0.1:tf;

%Linear Quadratic Controller - LQG --------------------------------------

%Creation of matrix A in xdot(t)=Ax(t)

A=zeros(2*N);
A=[zeros(N) eye(N); -inv(M)*K -inv(M)*D]; %Goes into LQG - system dynamics

C=3*[eye(N) zeros(N); zeros(N) zeros(N)];%obervability matrix %setting C=I reverts to LQR problem
Q=15*(C'*C); %Control Matrix - assumes equal control weight across system
B=zeros(2*N,1);
B(N+1:2*N)=inv(M)*Bbar;
R=1;
X=eye(2*N); %initial guess for ARE

pie_lqg=newt(A,B,Q,inv(R),X);  %steady state ARE solution
fbK_lqg=inv(R)*B'*pie_lqg;           %feedback K gain matrix 
pea_lqg = are(A',C'*C,B*B');
F_lqg = pea_lqg*C';
Ac_lqg = A-B*fbK_lqg-F_lqg*C;
A_lqg = [A  -B*fbK_lqg;F_lqg*C  Ac_lqg];

%disturbance rejection... but we don't have any right now.
%controller now is u(t) = -K*Xhat(t) + ufw(t)

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
        Sol(i,:)=Sol(i,:)+y(i,j).*ClampedBasis(j,:);
    end
end


mesh(xarr,t,Sol) %Use in presentation
axis([0 1 t0 tf -1 1])
title('Clamped Dynamic Beam With LQG Controller')
xlabel('length - x')
ylabel('time - t')
zlabel('Deflection - w(t,x)')
