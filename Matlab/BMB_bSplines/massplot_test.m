%Test script massplot

a=0;
b=2;
mass_length=0.25;
x_spacing=0.001;

x_temp=a:x_spacing:b;
t_temp=0:0.01:1;
f = @(x,t) sin(x);
for i=1:length(t_temp)
sol_matrix(i,:)=f(x_temp);
end

[x_mass,sol_mass,mp,xgap]=massplot(a,b,x_spacing,mass_length,sol_matrix);

mesh(x_mass(1:mp),t_temp,sol_mass(:,1:mp))
hold on
mesh(x_mass(mp+1:mp+xgap-1),t_temp,sol_mass(:,mp+1:mp+xgap-1),'EdgeColor','black')
mesh(x_mass(mp+xgap:length(x_mass)),t_temp,sol_mass(:,mp+xgap:length(x_mass)))