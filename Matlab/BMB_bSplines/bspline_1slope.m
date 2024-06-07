function [TestBasis,Basis,D_Basis,DD_Basis,ClampedBasis,D_ClampedBasis,DD_ClampedBasis,ClampedBMBBasis,D_ClampedBMBBasis,DD_ClampedBMBBasis,BMBBasis,D_BMBBasis,DD_BMBBasis,xarr,x_spacing] = cubicbspline2(a,b,num_splines)


h=(b-a)/(num_splines-1);
xi=0;
x1=linspace(-2*h,xi-h,6001);
x2=linspace(xi-h,xi,6001);
x3=linspace(xi,xi+h,6001);
x4=linspace(xi+h,xi+2*h,6000);
x=[x1 x2 x3 x4];

%Spline Indexing Variables
dxi=(b-a)/(num_splines-1); %spacing for xi
splines_begin=a-3*dxi;
splines_end=b+3*dxi;

%basic spline
y1 = @(x) 1/h^3*(x-(xi-2*h)).^3; 
y2 = @(x) 1/h^3.*(h^3+3*h^2.*(x-(xi-h))+3*h.*(x-(xi-h)).^2-3*(x-(xi-h)).^3);
y3 = @(x) 1/h^3.*(h^3+3*h^2.*((xi+h)-x)+3*h.*((xi+h)-x).^2-3*((xi+h)-x).^3);
y4 = @(x) 1/h^3.*((2*h+xi)-x).^3; 

%first derivative
y1prime= @(x) 1/h^3*(3*(2*h + x - xi).^2);
y2prime= @(x) 1/h^3*-3*(x - xi).*(4*h + 3*x - 3*xi);
y3prime= @(x) 1/h^3*(-3*(x - xi).*(4*h - 3*x + 3*xi));
y4prime= @(x) 1/h^3*(-3*(2*h - x + xi).^2);

%second derivative
y1pp= @(x) (6*(2*h + x - xi))*1/h^3;
y2pp= @(x) 1/h^3*-6*(2*h + 3*x - 3*xi);
y3pp= @(x) (-6*(2*h - 3*x + 3*xi))*1/h^3;
y4pp= @(x) (6*(2*h - x + xi))*1/h^3;


%y1first = @(x) y3-2*y2-2*y4;  %fix over whole spline
%dy1first = @(x) y3prime-2*y2prime-2*y4prime;
%ddy1first= @(x) y3pp-2*yy2pp -2*y4pp;
%yfirst=[y1first(x1)];
%dyfirst=[dy1first(x1)];
%ddyfirst=[ddy1first(x1)];

y = [y1(x1(1:length(x1)-1)) y2(x2(1:length(x2)-1)) y3(x3(1:length(x3)-1)) y4(x4)];
dy=[y1prime(x1(1:length(x1)-1)) y2prime(x2(1:length(x2)-1)) y3prime(x3(1:length(x3)-1)) y4prime(x4)];
ddy=[y1pp(x1(1:length(x1)-1)) y2pp(x2(1:length(x2)-1)) y3pp(x3(1:length(x3)-1)) y4pp(x4)];

%spline indexing

spline_length=length(y); %length of each spine - could be changed to offer more precision (also change lines 6-9)
interval_length=((num_splines-1)*spline_length+6*length(y))/4;  %calculating the length of the interval for spline storage
Bspline=zeros(num_splines+2,interval_length);    %initialize spline storage matrix
x_spacing=(splines_end-splines_begin)/(((num_splines-1)*spline_length+6*length(y))/4);
xarr=a:x_spacing:b; %array used for plotting splines

%Creating Bspline matrix (full splines N+2)
for i=1:spline_length
    for j=1:(num_splines+2)
        Bspline(j,i+(spline_length/4)*(j-1))=y(i);
    end
end

%Bspline basis vector matrix for [a,b]
Basis=zeros(num_splines,length(y)/4*(num_splines-1)); %Bbasis is Bspline with the two extra splines deleted
    
for i=3*length(y)/4+1:(length(y)/4*(num_splines-1)+3*length(y)/4+1)
    for j=1:(num_splines+2)
        Basis(j,i-3*length(y)/4)=Bspline(j,i);
    end
end

%D_BSpline - 1st derivative of spline functions (N+2)
D_Spline=zeros(13,interval_length);
for i=1:spline_length
    for j=1:(num_splines+2)
        D_Spline(j,i+(spline_length/4)*(j-1))=dy(i);
    end
end

%D_Basis - 1st derivative of basis functions for [a,b]
D_Basis=zeros(13,length(y)/4*(num_splines-1));
for i=3*length(y)/4+1:(length(y)/4*(num_splines-1)+3*length(y)/4+1)
    for j=1:(num_splines+2)
        D_Basis(j,i-3*length(y)/4)=D_Spline(j,i);
    end
end

%DD_Spline - 2nd derivative of spline functions (N+2)
DD_Spline=zeros(13,interval_length);
for i=1:spline_length
    for j=1:(num_splines+2)
        DD_Spline(j,i+(spline_length/4)*(j-1))=ddy(i);
    end
end

%DD_Basis - 2nd derivative of basis functions
DD_Basis=zeros(13,length(y)/4*(num_splines-1));
for i=3*length(y)/4+1:(length(y)/4*(num_splines-1)+3*length(y)/4+1)
    for j=1:(num_splines+2)
        DD_Basis(j,i-3*length(y)/4)=DD_Spline(j,i);
    end
end


%% Linear Combinations of Splines

%Clamped Beam left end condition (clamped at y=0, dy=0)
%B1=B2-2B1-2B3

%Clamped Basis fxn
ClampedBasis=zeros(num_splines,length(xarr));
ClampedBasis(1,:)=Basis(2,:)-2*Basis(1,:)-2*Basis(3,:); %linear combination
for i=4:(num_splines+2)
    ClampedBasis(i-2,:)=Basis(i,:);
end
%D_ClampedBasis - 1st derivative of clamped basis fxn
D_ClampedBasis=zeros(num_splines,length(xarr));
D_ClampedBasis(1,:)=D_Basis(2,:)-2*D_Basis(1,:)-2*D_Basis(3,:); %linear combination
for i=4:(num_splines+2)
    D_ClampedBasis(i-2,:)=D_Basis(i,:);
end
%DD_ClampedBasis - 2nd derivative of clamped basis fxn
DD_ClampedBasis=zeros(num_splines,length(xarr));
DD_ClampedBasis(1,:)=DD_Basis(2,:)-2*DD_Basis(1,:)-2*DD_Basis(3,:); %linear combination
for i=4:(num_splines+2)
    DD_ClampedBasis(i-2,:)=DD_Basis(i,:);
end


%Clamped Left Beam with a zero slope condition at the mass location w/
%free end condition
ClampedBMBBasis=zeros(num_splines,length(xarr));
ClampedBMBBasis((num_splines-1)/2,:)=ClampedBasis((num_splines-1)/2,:)-2*ClampedBasis((num_splines-1)/2-1,:)-2*ClampedBasis((num_splines-1)/2+1,:);
for i=1:((num_splines-1)/2-1)
    ClampedBMBBasis(i,:)=ClampedBasis(i,:);
end
for i=((num_splines-1)/2+1):num_splines
    ClampedBMBBasis(i,:)=ClampedBasis(i,:);
end

D_ClampedBMBBasis=zeros(num_splines,length(xarr));
D_ClampedBMBBasis((num_splines-1)/2,:)=D_ClampedBasis((num_splines-1)/2,:)-2*D_ClampedBasis((num_splines-1)/2-1,:)-2*D_ClampedBasis((num_splines-1)/2+1,:);
for i=1:((num_splines-1)/2-1)
    D_ClampedBMBBasis(i,:)=D_ClampedBasis(i,:);
end
for i=((num_splines-1)/2+1):num_splines
    D_ClampedBMBBasis(i,:)=D_ClampedBasis(i,:);
end


DD_ClampedBMBBasis=zeros(num_splines,length(xarr));
DD_ClampedBMBBasis((num_splines-1)/2,:)=DD_ClampedBasis((num_splines-1)/2,:)+DD_ClampedBasis((num_splines-1)/2-1,:)+DD_ClampedBasis((num_splines-1)/2+1,:);
for i=1:((num_splines-1)/2-1)
    DD_ClampedBMBBasis(i,:)=DD_ClampedBasis(i,:);
end
for i=((num_splines-1)/2+1):num_splines
    DD_ClampedBMBBasis(i,:)=DD_ClampedBasis(i,:);
end
ClampedBMBBasis(((num_splines-1)/2-1),:)=[];
ClampedBMBBasis(((num_splines-1)/2),:)=[];
D_ClampedBMBBasis(((num_splines-1)/2-1),:)=[];
D_ClampedBMBBasis(((num_splines-1)/2),:)=[];
DD_ClampedBMBBasis(((num_splines-1)/2-1),:)=[];
DD_ClampedBMBBasis(((num_splines-1)/2),:)=[];



% BMB Basis Vectors with zero slope condition. Unclamped.
BMBBasis=Basis;
BMBBasis((num_splines+3)/2,:)=2/3*(Basis((num_splines+3)/2+1,:)+Basis((num_splines+3)/2-1,:)+Basis((num_splines+3)/2,:));
BMBBasis((num_splines+3)/2-1,:)=[];
BMBBasis((num_splines+3)/2,:)=[];

D_BMBBasis=D_Basis;
D_BMBBasis((num_splines+3)/2,:)=D_Basis((num_splines+3)/2+1,:)+D_Basis((num_splines+3)/2-1,:)+D_Basis((num_splines+3)/2,:);
D_BMBBasis((num_splines+3)/2-1,:)=[];
D_BMBBasis((num_splines+3)/2,:)=[];

DD_BMBBasis=DD_Basis;
DD_BMBBasis((num_splines+3)/2,:)=DD_Basis((num_splines+3)/2+1,:)+DD_Basis((num_splines+3)/2-1,:)+DD_Basis((num_splines+3)/2,:);
DD_BMBBasis((num_splines+3)/2-1,:)=[];
DD_BMBBasis((num_splines+3)/2,:)=[];

%Uncomment to force slope of 1 at mass location
TestBasis=BMBBasis;
 D_TestBasis=BMBBasis;
 DD_TestBasis=BMBBasis;
 %TestBasis((num_splines+1)/2,1:(length(xarr)-1)/2+1)=1/12*Basis((num_splines+3)/2-1,1:(length(xarr)-1)/2+1)-2/21*Basis((num_splines+3)/2,1:(length(xarr)-1)/2+1);%+Basis((num_splines+3)/2,:);
 TestBasis((num_splines+1)/2,:)=-13/12*Basis((num_splines+3)/2-1,:)+Basis((num_splines+3)/2+1,:);
 %BMBBasis((num_splines+1)/2+1,:)=h/24*BMBBasis((num_splines+1)/2-1,:)-h/24*BMBBasis((num_splines+1)/2+1,:);

 %D_BMBBasis((num_splines+1)/2-1,:)=-D_BMBBasis((num_splines+1)/2-1,:);%-h/24*D_BMBBasis((num_splines+1)/2+1,:);
 %D_BMBBasis((num_splines+1)/2+1,:)=h/24*D_BMBBasis((num_splines+1)/2-1,:)-h/24*D_BMBBasis((num_splines+1)/2+1,:);

 %DD_BMBBasis((num_splines+3)/2,:)=-DD_BMBBasis((num_splines+1)/2-1,:);%-h/24*DD_BMBBasis((num_splines+1)/2+1,:);
 %DD_BMBBasis((num_splines+1)/2+1,:)=h/24*DD_BMBBasis((num_splines+1)/2-1,:)-h/24*DD_BMBBasis((num_splines+1)/2+1,:);








