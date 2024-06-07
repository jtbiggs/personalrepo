function [Xprime] = dbeam_lqg(t,X,A,B,K,ufw,F,d,Xend)
% u =-Kx
%Xdot = Ax + Bu + Fd
%need u(t) = -K*xhat(t) + ufw(t)  since we're using the state estimate
Xprime = A*X-A*Xend;
end