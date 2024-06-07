%Redo of the time invariant BMB model with clamped left condition over the
%interval [0,2]

%IC: y(0,x)=-0.5x
%Forcing function: y''''=2pi^2sin(pix)

%Static Variables
a=0;    %beginning of interval
b=2;    %end of interval
f = @(x) 2*pi^2*sin(pi*x); %forcing function
N=23; %N-2 total splines for the interval [a,b] or (N-3)/2 splines per beam with one shared middle spline
mass_length=0.25;
% K1=zeros(((N-3)/2-1));
% K2=zeros((N-3)/2);


%import basis functions
[Basis,D_Basis,DD_Basis,ClampedBasis,D_ClampedBasis,DD_ClampedBasis,ClampedBMBBasis,D_ClampedBMBBasis,DD_ClampedBMBBasis,xarr,x_spacing] = cubicbspline(a,b,N);
N=N-2;
K=zeros(N);
for i=1:N
     for j=1:N
     K(i,j)=trapz(xarr,DD_ClampedBMBBasis(i,:).*DD_ClampedBMBBasis(j,:));
     end
end

F=zeros(N,1);
for i=1:N
     F(i)=trapz(xarr,f(xarr).*ClampedBMBBasis(i,:));
end
C=zeros(N,1);
C=K\F;
Sol=zeros(1,length(xarr));

for i=1:N
    Sol=Sol+C(i).*ClampedBMBBasis(i,:);
end
[xmass,Solmass,mp,xgap] = massplot(a,b,x_spacing,mass_length,Sol);

plot(xmass(1:mp),Solmass(1:mp),'LineWidth',2)
hold on
plot(xmass(mp+1:mp+xgap-1),Solmass((mp+1:mp+xgap-1)),'color','#b4b4b4','LineWidth',6)
plot(xmass(mp+xgap:length(xmass)),Solmass(mp+xgap:length(xmass)),'LineWidth',2)
grid
axis([0 2.25 -2 2])
xlabel('Length - x')
ylabel('Deflection - w(x)')
title('ODE BMB System')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %K1 stiffness matrix (left beam)
% for i=1:((N-3)/2-1)
%     for j=1:((N-3)/2-1)
%     K1(i,j)=trapz(xarr,DD_ClampedBMBBasis(i,:).*DD_ClampedBMBBasis(j,:));
%     end
% end
% 
% %K2 stiffness matrix
% for i=((N-3)/2+3):N
%     for j=((N-3)/2+3):N
%         K2(i-((N-3)/2+2),j-((N-3)/2+2))=trapz(xarr,DD_ClampedBMBBasis(i,:).*DD_ClampedBMBBasis(j,:));
%     end
% end

% %K matrix - combination of K1 & K2
% K=zeros((N-6));
% K1_zeros=zeros((N-6));
% K2_zeros=zeros((N-6));
% 
% for i=1:((N-3)/2-1)
%     for j=1:((N-3)/2-1)
%         K1_zeros(i,j)=K1(i,j);
%     end
% end
% for i=1:((N-3)/2+1)
%     for j=1:((N-3)/2+1)
%         K2_zeros(i+(N-3)/2-4,j+(N-3)/2-4)=K2(i,j);
%     end
% end
% K=K1_zeros+K2_zeros;
% 
% %Integral of forcing function over left and right beams
% F=zeros(1,N-6);
% F1=zeros(N,1);
% F2=zeros((N-3)/2+1,1);
% F1_zeros=zeros(N-6,1);
% F2_zeros=zeros(N-6,1);
% 
% for i=1:((N-3)/2-1)
%     F1(i)=trapz(xarr,f(xarr).*ClampedBMBBasis(i,:));
%     F1_zeros(i)=F1(i);
% end
% for i=((N-3)/2+3):N
%     F2(i-((N-3)/2+2))=trapz(xarr,f(xarr).*ClampedBMBBasis(i,:));
%     F2_zeros(i-((N-3)/2-4))=F2(i-((N-3)/2+2));
% end
% F=F1_zeros+F2_zeros;
% 
% %Coefficient Vector
% C=zeros(N-6,1);
% C=inv(K)*F;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
