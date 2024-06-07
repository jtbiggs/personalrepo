function [x_mass,sol_mass,mp,xgap] = massplot(a,b,x_spacing,mass_length,sol_matrix,num_splines,t)
%Generates fixed y-values at the mass location for BMB solutions
N=num_splines;
b_mass=b+mass_length;
x=a:x_spacing:b;
x_mass=a:x_spacing:b_mass;
t_length=length(sol_matrix(:,1));
x_length=length(x_mass);
sol_mass=zeros(t_length,x_length);
Solplot=zeros(t_length,N);
xgap=length(b:x_spacing:b_mass);

mp=(length(x)-1)/2+1;
y=sol_matrix;
for i=1:t_length
sol_mass(i,1:mp)=y(i,1:mp);
sol_mass(i,(mp+1):(mp+xgap))=y(i,mp);
sol_mass(i,mp+xgap:length(x_mass))=y(i,mp+1:length(x));
Solplot(i,:)=sol_mass(i,1:ceil((length(x_mass)+1)/N):end);
end
Solplot=[Solplot sol_mass(:,end)];
xplot=x_mass(1:ceil((length(x_mass)+1)/N):end);
xplot=[xplot x_mass(end)];

mesh(xplot,t,Solplot)
hold on
mesh(x_mass(mp+1:mp+xgap-1),t,sol_mass(:,(mp+1:mp+xgap-1)),'EdgeColor','#b4b4b4')
grid on
