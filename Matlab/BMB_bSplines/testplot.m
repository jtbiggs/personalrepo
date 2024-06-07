xplot=zeros(length(t),9);
Solplot=zeros(length(t),9);
for i=1:length(t)
Solplot(i,:)=Solmass(i,1:(length(xmass)+1)/N:end);
end
Solplot=[Solplot Solmass(:,end)];
xplot=xmass(1:(length(xmass)+1)/N:end);
xplot=[xplot xmass(end)];
mesh(xplot,t,Solplot)
x2=xmass(mp:)