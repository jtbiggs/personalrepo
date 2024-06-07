function findb = findb(t,b,flag,A,B,invR,P,track_state)
    t;
        findb = -(A-B*invR*B'*P)'*b; % - P*A*track_state;