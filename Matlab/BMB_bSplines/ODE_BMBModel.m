%ODE BMB Model
%% Initialization
%constants & initial conditions
a=0;
b=1;
a1=1;
b2=2;
n=10;
slength=4000;
nsl=n*slength/4;
nsl2=(n+3)*slength/4; %(n+3)
f1 = @(x) 2*pi*sin(pi*x); 
f2 = f1;
F1=zeros(2*(n-1)+1,1);
F2=zeros(2*(n-1)+1,1);
c=zeros(2*(n),1);
Sol=0;
%spline functions
dx_spline=(b-a)/(nsl-1);
dx_spline2=(b2-a1)/(nsl2-1);
xarr1=a:dx_spline:b;
x3=linspace(0,1,nsl);
xarr2=1:dx_spline:2; %FIX 20,000 to 23.000
xarr2_plot=1:dx_spline2:2;
%initialize storage matrices
Smat=zeros(n-1,nsl);
Smat2=zeros(n,nsl);
Ds=zeros(n-1,nsl);
DDs=zeros(n-1,nsl);
Ds2=zeros(n,nsl);
DDs2=zeros(n,nsl);
K1=zeros(n-1);%empty stiffness matrix 1
K2=zeros(n+1);%empty stiffness matrix 2
M=zeros(2*(n+1)); %combined matrix
%endpoint spacing for splines
dx=(b-a)/n;
q=a-3*dx;
p=a+dx;

%% Generation of Cubic B-Spline Basis Functions

y=zeros(n+3,slength);
dy=zeros(n+3,slength);

%linear combination spline clamped left beam
weights=[1 -2 -2];
[s1,X,s1dy,s1ddy]=cubicbspline(a-3*dx,a+dx);
lincombo=weights(:)*cubicbspline(a-3*dx,a+dx);
dylincombo=weights(:)*s1dy;
ddylincombo=weights(:)*s1ddy;
x1=linspace(0,3*dx,slength-(slength/4));
totcombo=zeros(slength-(slength/4));

%1st Spline setup
com1=lincombo(2,0.75*slength+1:slength)+lincombo(1,.5*slength+1:0.75*slength)+lincombo(3,0.25*slength+1:0.5*slength);
com2=lincombo(1,0.75*slength+1:slength)+lincombo(3,.5*slength+1:0.75*slength);
com3=lincombo(3,0.75*slength+1:slength);
totcombo=[com1 com2 com3];

%1st derivative of 1st spline setup
dcom1=dylincombo(2,0.75*slength+1:slength)+dylincombo(1,.5*slength+1:0.75*slength)+dylincombo(3,0.25*slength+1:0.5*slength);
dcom2=dylincombo(1,0.75*slength+1:slength)+dylincombo(3,.5*slength+1:0.75*slength);
dcom3=dylincombo(3,0.75*slength+1:slength);
dtotcombo=[dcom1 dcom2 dcom3];

%2nd derivative of 1st spline setup
ddcom1=ddylincombo(2,0.75*slength+1:slength)+ddylincombo(1,.5*slength+1:0.75*slength)+ddylincombo(3,0.25*slength+1:0.5*slength);
ddcom2=ddylincombo(1,0.75*slength+1:slength)+ddylincombo(3,.5*slength+1:0.75*slength);
ddcom3=ddylincombo(3,0.75*slength+1:slength);
ddtotcombo=[ddcom1 ddcom2 ddcom3];

%Linear combination y'_L(b)=y'_R(b)
%Zero Derivative at b with continuous position
h=dx;
xi1=b-h;
%Spline 1
Sp1_1 = @(x) 1/h^3*(x-(xi1-2*h)).^3;
Sp1_2 = @(x) 1/h^3.*(h^3+3*h^2.*(x-(xi1-h))+3*h.*(x-(xi1-h)).^2-3*(x-(xi1-h)).^3);
Sp1_3 = @(x) 1/h^3.*(h^3+3*h^2.*((xi1+h)-x)+3*h.*((xi1+h)-x).^2-3*((xi1+h)-x).^3);
Sp1_4 = @(x) 1/h^3.*((2*h+xi1)-x).^3;

DSp1_1= @(x) 1/h^3*(3*(2*h + x - xi1).^2);
DSp1_2= @(x) 1/h^3*-3*(x - xi1).*(4*h + 3*x - 3*xi1);
DSp1_3= @(x) 1/h^3*(-3*(x - xi1).*(4*h - 3*x + 3*xi1));
DSp1_4= @(x) 1/h^3*(-3*(2*h - x + xi1).^2);

DDSp1_1= @(x) (6*(2*h + x - xi1))*1/h^3;
DDSp1_2= @(x) 1/h^3*-6*(2*h + 3*x - 3*xi1);
DDSp1_3= @(x) (-6*(2*h - 3*x + 3*xi1))*1/h^3;
DDSp1_4= @(x) (6*(2*h - x + xi1))*1/h^3;
%Spline 2
xi2=b;
Sp2_1 = @(x) 1/h^3*(x-(xi2-2*h)).^3;
Sp2_2 = @(x) 1/h^3.*(h^3+3*h^2.*(x-(xi2-h))+3*h.*(x-(xi2-h)).^2-3*(x-(xi2-h)).^3);
Sp2_3 = @(x) 1/h^3.*(h^3+3*h^2.*((xi2+h)-x)+3*h.*((xi2+h)-x).^2-3*((xi2+h)-x).^3);
Sp2_4 = @(x) 1/h^3.*((2*h+xi2)-x).^3;

DSp2_1= @(x) 1/h^3*(3*(2*h + x - xi2).^2);
DSp2_2= @(x) 1/h^3*-3*(x - xi2).*(4*h + 3*x - 3*xi2);
DSp2_3= @(x) 1/h^3*(-3*(x - xi2).*(4*h - 3*x + 3*xi2));
DSp2_4= @(x) 1/h^3*(-3*(2*h - x + xi2).^2);

DDSp2_1= @(x) (6*(2*h + x - xi2))*1/h^3;
DDSp2_2= @(x) 1/h^3*-6*(2*h + 3*x - 3*xi2);
DDSp2_3= @(x) (-6*(2*h - 3*x + 3*xi2))*1/h^3;
DDSp2_4= @(x) (6*(2*h - x + xi2))*1/h^3;

%Spline 3
xi3=b+h;
Sp3_1= @(x) 1/h^3*(x-(xi3-2*h)).^3;
Sp3_2 = @(x) 1/h^3.*(h^3+3*h^2.*(x-(xi3-h))+3*h.*(x-(xi3-h)).^2-3*(x-(xi3-h)).^3);
Sp3_3 = @(x) 1/h^3.*(h^3+3*h^2.*((xi3+h)-x)+3*h.*((xi3+h)-x).^2-3*((xi3+h)-x).^3);
Sp3_4 = @(x) 1/h^3.*((2*h+xi3)-x).^3;

DSp3_1= @(x) 1/h^3*(3*(2*h + x - xi3).^2);
DSp3_2= @(x) 1/h^3*-3*(x - xi3).*(4*h + 3*x - 3*xi3);
DSp3_3= @(x) 1/h^3*(-3*(x - xi3).*(4*h - 3*x + 3*xi3));
DSp3_4= @(x) 1/h^3*(-3*(2*h - x + xi3).^2);

DDSp3_1= @(x) (6*(2*h + x - xi3))*1/h^3;
DDSp3_2= @(x) 1/h^3*-6*(2*h + 3*x - 3*xi3);
DDSp3_3= @(x) (-6*(2*h - 3*x + 3*xi3))*1/h^3;
DDSp3_4= @(x) (6*(2*h - x + xi3))*1/h^3;

%Linear Combination at mass location
totcombo2=zeros(1,slength+slength/2);
dtotcombo2=zeros(1,slength+slength/2);
ddtotcombo2=zeros(1,slength+slength/2);

com2_1=[Sp1_1(linspace(b-3*h,b-2*h,slength/4))];
com2_2=[Sp1_2(linspace(b-2*h,b-h,slength/4))+Sp2_1(linspace(b-2*h,b-h,slength/4))];
com2_3=[Sp1_3(linspace(b-h,b,slength/4))+Sp3_1(linspace(b-h,b,slength/4))+Sp2_2(linspace(b-h,b,slength/4))];
com2_4=[Sp1_4(linspace(b,b+h,slength/4))+Sp3_2(linspace(b,b+h,slength/4))+Sp2_3(linspace(b,b+h,slength/4))];
com2_5=[Sp3_3(linspace(b+h,b+2*h,slength/4))+Sp2_4(linspace(b+h,b+2*h,slength/4))];
com2_6=[Sp3_4(linspace(b+2*h,b+3*h,slength/4))];

Dcom2_1=[DSp1_1(linspace(b-3*h,b-2*h,slength/4))];
Dcom2_2=[DSp1_2(linspace(b-2*h,b-h,slength/4))+DSp2_1(linspace(b-2*h,b-h,slength/4))];
Dcom2_3=[DSp1_3(linspace(b-h,b,slength/4))+DSp3_1(linspace(b-h,b,slength/4))+DSp2_2(linspace(b-h,b,slength/4))];
Dcom2_4=[DSp1_4(linspace(b,b+h,slength/4))+DSp3_2(linspace(b,b+h,slength/4))+DSp2_3(linspace(b,b+h,slength/4))];
Dcom2_5=[DSp3_3(linspace(b+h,b+2*h,slength/4))+DSp2_4(linspace(b+h,b+2*h,slength/4))];
Dcom2_6=[DSp3_4(linspace(b+2*h,b+3*h,slength/4))];

DDcom2_1=[DDSp1_1(linspace(b-3*h,b-2*h,slength/4))];
DDcom2_2=[DDSp1_2(linspace(b-2*h,b-h,slength/4))+DDSp2_1(linspace(b-2*h,b-h,slength/4))];
DDcom2_3=[DDSp1_3(linspace(b-h,b,slength/4))+DDSp3_1(linspace(b-h,b,slength/4))+DDSp2_2(linspace(b-h,b,slength/4))];
DDcom2_4=[DDSp1_4(linspace(b,b+h,slength/4))+DDSp3_2(linspace(b,b+h,slength/4))+DDSp2_3(linspace(b,b+h,slength/4))];
DDcom2_5=[DDSp3_3(linspace(b+h,b+2*h,slength/4))+DDSp2_4(linspace(b+h,b+2*h,slength/4))];
DDcom2_6=[DDSp3_4(linspace(b+2*h,b+3*h,slength/4))];


totcombo2=[com2_1 com2_2 com2_3 com2_4 com2_5 com2_6];
dtotcombo2=[Dcom2_1 Dcom2_2 Dcom2_3 Dcom2_4 Dcom2_5 Dcom2_6];
ddtotcombo2=[DDcom2_1 DDcom2_2 DDcom2_3 DDcom2_4 DDcom2_5 DDcom2_6];

%Generation of remaining left splines
for k=1:n+3
    [y(k,:),x,dy(k,:),ddy(k,:)]=cubicbspline(q,p);
    q=q+dx;
    p=p+dx;
end
%Generation of right splines
r=a1-3.*dx;
s=a1+dx;
for k=1:n+3
    [y2(k,:),x,dy2(k,:),ddy2(k,:)]=cubicbspline(r,s);
    r=r+dx;
    s=s+dx;
end

%storage of splines in matrix 1
for i=1:(slength-slength/4)
    Smat(1,i)=totcombo(1,i); 
    Ds(1,i)=dtotcombo(1,i);
    DDs(1,i)=ddtotcombo(1,i);
end
for j=1:slength
    Smat(2,j)=y(1,j);
    Ds(2,j)=dy(1,j);
    DDs(2,j)=ddy(1,j);
end
for k=2:n-3
    for p=1:slength
    Smat(k+1,p+slength/4*(k-1))=y(k,p);
    Ds(k+1,p+slength/4*(k-1))=dy(k,p);
    DDs(k+1,p+slength/4*(k-1))=ddy(k,p);
    end
end
for i=1:(6*slength/4)
    Smat(n-1,i+((n-3)*(slength/4)-1))=totcombo2(1,i);
    Ds(n-1,i+((n-3)*(slength/4)-1))=dtotcombo2(1,i);
    DDs(n-1,i+((n-3)*(slength/4)-1))=ddtotcombo2(1,i);
end


for i=0.5*length(totcombo2):(slength+slength/2)
    Smat2(1,i-(0.5*length(totcombo2)-1))=totcombo2(1,i);
    Ds2(1,i-(0.5*length(totcombo2)-1))=dtotcombo2(1,i);
    DDs2(1,i-(0.5*length(totcombo2)-1))=ddtotcombo2(1,i);
end
% for i=slength/4:slength
%     Smat2(2,i-slength/4+1)=y(1,i);
%     Ds2(2,i-slength/4+1)=dy(1,i);
%     DDs2(2,i-slength/4+1)=ddy(1,i);
% end

for i=2:n+1
    for j=1:slength
        Smat2(i,j+slength/4*(i-2))=y(i,j);
        Ds2(i,j+slength/4*(i-2))=dy(i,j);
        DDs2(i,j+slength/4*(i-2))=ddy(i,j);
    end
end

Smat=Smat(:,1:nsl);
Ds=Ds(:,1:nsl);
DDs=DDs(:,1:nsl);
Smat2=Smat2(:,1:nsl);
Ds2=Ds2(:,1:nsl);
Ds2=Ds2(:,1:nsl);

%Finite Element Method (ODE)
 for i=1:n-1
     for j=1:n-1
         K1(i,j)=trapz(dx_spline,DDs(i,:).*DDs(j,:));
     end
 end
 for i=1:n+1
     for j=1:n+1
         K2(i,j)=trapz(dx_spline,DDs2(i,:).*DDs2(j,:));
     end
 end 
M1=[zeros(2*n-1)];
M2=[zeros(2*n-1)];
M=[zeros(2*n-1)];
for i=1:length(K1(:,1))
    for j=1:length(K1(:,1))
    M1(i,j)=K1(i,j);
    end
end
for i=1:length(K2(:,1))
    for j=1:length(K2(:,1))
        M2(i+(n-2),j+(n-2))=K2(i,j);
    end
end
M=M1+M2;
for i=1:n-1
F1(i)=trapz(xarr1,f1(xarr1).*Smat(i,:));
end
for i=1:n+1
F2(n+i-2)=trapz(xarr2, f2(xarr2).*Smat2(i,:));
end
F=zeros(n-1,1);
F=F1+F2;
c=inv(M)*F;
Sol1=0;
Sol2=Sol1;
for i=1:n-1
Sol1=Sol1+c(i)*Smat(i,:);
end
for i=n:2*n-1
Sol2=Sol2+c(i)*Smat2(i-(n-1),:);
end
x4=linspace(0.7,2,13000);

figure(1)
% plot(0:dx_spline:1,Sol1)
% grid on
% hold on
% plot(1:dx_spline:2,Sol2)
%hold on
%plot(xarr1,Smat)
% hold on
%figure(2)
%plot(xarr2,Smat2)

plot(xarr1,DDs)
% hold on
figure(2)
plot(xarr2_plot,DDs2)