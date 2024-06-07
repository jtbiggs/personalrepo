function xprime = dbeam_lqr(t,x,A,K,B)
    %u=-Kx 
    xprime= A*x+B*(-K*x); %xdot=Ax