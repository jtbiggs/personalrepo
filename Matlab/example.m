a=12;
b=45;
A=[];
for i=1:14
    A=[A i];
end
x=linspace(0,12,100);
y= @(x) 0.5*x-4;

fplot3(@(t) cos(t),@(t) sin(t),@(t) t)
grid on