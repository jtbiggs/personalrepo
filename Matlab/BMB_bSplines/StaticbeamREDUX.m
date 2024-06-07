%Redo of time-independent clamped beam model (ODE)
clc

%Fixed Variables
a=0;
b=2;
N=11;
K=zeros(N);
F=zeros(N,1);
C=zeros(N,1);

f = @(x) 2*pi^2*sin(pi*x); %forcing function on beam

%Import cubic b-splines from function
[Basis,D_Basis,DD_Basis,ClampedBasis,D_ClampedBasis,DD_ClampedBasis,~,~,~,xarr]=cubicbspline(a,b,N);
%Only need Basis and DD_Spline for this example

%Generation of K stiffness matrix for finite element method

for i=1:N
    for j=1:N
    K(i,j)=simps(xarr,DD_ClampedBasis(i,:).*DD_ClampedBasis(j,:));
    end
end

%Integral of forcing function over [a,b]
for i=1:N
    F(i,:)=trapz(xarr,f(xarr).*ClampedBasis(i,:));
end
Kreal=[81000	12000	750	-1500	0	0	0	0	0	0	0;
       12000	12000	-6750	0	750	0	0	0	0	0	0;
       750	-6750	12000	-6750	0	750	0	0	0	0	0;
      -1500	0	-6750	12000	-6750	0	750	0	0	0	0;
        0	750	0	-6750	12000	-6750	0	750	0	0	0;
        0	0	750	0	-6750	12000	-6750	0	750	0	0;
        0	0	0	750	0	-6750	12000	-6750	0	750	0;
        0	0	0	0	750	0	-6750	12000	-6750	0	750;
        0	0	0	0	0	750	0	-6750	10500	-4500	0;
        0	0	0	0	0	0	750	0	-4500	6000	-2250;
        0	0	0	0	0	0	0	750	0	-2250	1500];
%Generation of Solution coefficient vector
C=K\F;
Creal=K\F;
%ODE solution vector - sum of the product of sol coef. and clamped basis
Sol=zeros(1,length(xarr));
for i=1:N
    Sol=Sol+C(i).*ClampedBasis(i,:);
end
   
%Plot of solution

%Exact Solution
Sol_exact= (-6*pi*xarr - 6*pi^3*xarr.^2+pi^3*xarr.^3 + 6*sin(pi*xarr))/(3*pi^2);

plot(xarr,Sol,'--o','MarkerIndices',1:(length(Sol)+4)/N:length(Sol))
hold on
plot(xarr,Sol_exact)
axis([a b -20 20]);
title('ODE Clamped Beam Solutions')
xlabel('Length - x')
ylabel('Beam Deflection - y(x)')
legend({'Numerical Solution','Exact Solution'},'Location','northeast')




