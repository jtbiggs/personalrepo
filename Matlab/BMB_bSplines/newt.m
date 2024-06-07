%This file uses the Newton-Kleinman method to find a (hopefully) positive semidefinite solution to
%the ARE.  See Datta p. 567 for algorithm.
%Given an initial condition X, and matrices A, B, Q, invR, this solves the following ARE:
%A'P+PA+Q-PB(invR)B'P

function P = newt(A,B,Q,invR,X)
S = B*invR*B';

rel_chng = 1e-3;
Newt_Klein_conv_check = 1;
while Newt_Klein_conv_check>rel_chng    
    A_i = A - S*X;
    R_c = A'*X + X*A + Q - X*S*X;
    delta = lyap(A_i',R_c);
    X = X + delta;
    Newt_Klein_conv_check = norm(X-(X-delta))/norm(X-delta);
        if Newt_Klein_conv_check<=rel_chng
            soln=X;
        end 
end

P=soln;

