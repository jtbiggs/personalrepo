%Exact Solution BMB_ODE redux

%y''''=2*pi*sin(pi*x)
%condition: y(0)=y'(0)=y'(1)=y''(1)=y'''(1)=0

syms f(x) X Y
Df=diff(f,x);
Df2=diff(f,x,2);
Df3=diff(f,x,3);
Df4=diff(f,x,4);
ode = Df4==2*pi^2*sin(pi*x);
cond1 = f(0) == 0;
cond2 = Df(0) == 0;
cond3 = Df2()

%conds = [cond1 cond2];

sol=dsolve(ode);




% [VF,Sbs]=odeToVectorField(ode);
% odefxn=matlabFunction(VF,'Vars',{x,Y});
% 
% [X,Y]=ode45(odefxn,[0 1],[0 0 0 0]);
% figure
% plot(X,Y)
% grid
% legend(string(Sbs))