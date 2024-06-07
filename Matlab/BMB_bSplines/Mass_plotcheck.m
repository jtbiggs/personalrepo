
%plot check
a=0;
b=2;
b2=2.25;
h=0.001;
x=a:h:b;
x_adj=a:h:b2;
y_adj=zeros(1,length(x_adj));
xgap=length(b:h:b2);
f = @(x) x.^2; 
y=f(x);
midpoint=(length(x)-1)/2+1;
y_adj(1:midpoint)=y(1:midpoint);
y_adj((midpoint+1):(midpoint+xgap))=y(midpoint);
y_adj(midpoint+xgap:length(x_adj))=y(midpoint+1:length(x));
% plot(x_adj(1:midpoint),y_adj(1:midpoint))
% grid on
% hold on
% plot(x_adj((midpoint+1):(midpoint+xgap)),y_adj((midpoint+1):(midpoint+xgap)),'LineWidth',10)
% hold on
% plot(x_adj(midpoint+xgap:length(x_adj)),y_adj(midpoint+xgap:length(x_adj)))
plot(x_adj,y_adj)
axis([0 b2 0 4])