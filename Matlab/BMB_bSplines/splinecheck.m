
clc
a=0;
b=2;
N=11; %this generates N splines
[TestBasis,Basis,D_Basis,DD_Basis,ClampedBasis,D_ClampedBasis,DD_ClampedBasis,ClampedBMBBasis,D_ClampedBMBBasis,DD_ClampedBMBBasis,BMBBasis,D_BMBBasis,DD_BMBBasis,xarr,~]=bspline_1slope(a,b,N);
%[~,~,~,~,~,~,y,~,~,x,~]=cubicbspline(a,b,N);
figure
%x=linspace(4,8,24000);
plot(xarr,TestBasis)
title('f b-spline BMB Basis')
axis([0 2 -8 5]);
%xticks([4 5 6 7 8])
%xticklabels({'x_{i-2}','x_{i-1}','x_i','x_{i+1}','x_{i+2}'})
grid on
% figure
% plot(xarr,D_ClampedBasis)
% grid onasdf 

% Kreal=[81000	12000	750	-1500	0	0	0	0	0	0	0;
%        12000	12000	-6750	0	750	0	0	0	0	0	0;
%        750	-6750	12000	-6750	0	750	0	0	0	0	0;
%       -1500	0	-6750	12000	-6750	0	750	0	0	0	0;
%         0	750	0	-6750	12000	-6750	0	750	0	0	0;
%         0	0	750	0	-6750	12000	-6750	0	750	0	0;
%         0	0	0	750	0	-6750	12000	-6750	0	750	0;
%         0	0	0	0	750	0	-6750	12000	-6750	0	750;
%         0	0	0	0	0	750	0	-6750	10500	-4500	0;
%         0	0	0	0	0	0	750	0	-4500	6000	-2250;
%         0	0	0	0	0	0	0	750	0	-2250	1500];
