function [exact_sol] = DynamicbeamExact(a,b,N,t)

[~,~,~,~,~,~,~,~,~,xarr]=cubicbspline(a,b,N);

%Solve cosh(B)+1/cos(B)=0 for the 1st 10 solutions to get Beta_n
syms B
eqn = cosh(B)*cos(B)+1==0; 
beta1=double(vpasolve(eqn,B,1.5));
beta2=double(vpasolve(eqn,B,4.6));
beta3=double(vpasolve(eqn,B,7.855));
beta4=double(vpasolve(eqn,B,10.99));
beta5=double(vpasolve(eqn,B,[14.13 14.14]));
beta6=double(vpasolve(eqn,B,[17.27 17.28]));
beta7=double(vpasolve(eqn,B,[20.42 20.43]));
beta8=double(vpasolve(eqn,B,[23.56 23.57]));
beta9=double(vpasolve(eqn,B,[26.7 26.71]));
beta10=double(vpasolve(eqn,B,[29.84 29.85]));

Beta=[beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 beta9 beta10];

f = -0.5.*xarr;
X=zeros(10,length(xarr));
for i=1:10
    X(i,:)=cos(Beta(i).*xarr)-cosh(Beta(i).*xarr)-((cosh(Beta(i))+cos(Beta(i)))/(sinh(Beta(i))+sin(Beta(i))))*(sin(Beta(i).*xarr)-sinh(Beta(i).*xarr));
end
X(:,1)=zeros(10,1); 
kappa=zeros(10,1);
for i=1:10
kappa(i)=trapz(xarr,X(i,:).^2);
end

A=zeros(10,1);
for i=1:10
    A(i)=1./kappa(i).*trapz(xarr,f.*X(i,:));
end

exact_sol=zeros(length(t),length(xarr));
for i=1:length(t)
    for j=1:10
    exact_sol(i,:)=exact_sol(i,:)+X(j,:).*A(j).*cos(Beta(j).^2.*t(i));
    end
end
