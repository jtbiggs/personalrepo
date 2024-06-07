function z = liftfun(t,x,BMBBasis,xarr,A,G,k1,k2,k3,k4,k5,liftc,u,N,M)
    %This function builds Fbar matrix
  
    Fbar=zeros(N,1);

    %Probably not needed.... probably
    %Remove middle spline from basis since point mass has no lift
    %LiftBasis=BMBBasis;
    %LiftBasis((num_bsplines+1)/2,:)=zeros(1,length(1,BMBBasis));
    disp('t= ')
    disp(t)
    %integral for system Fbar matrix (left and right beams)
    for i=1:N
        Fbar(i)=simps(xarr,liftc*(k1+k2*sin(k3*atan((x(N+1:2*N,1)'*BMBBasis+k5)/u)-k4)).*BMBBasis(i,:));
    end
    
    F=zeros(2*N,1);
    F(N+1:2*N,1)=M\Fbar;
    
    z=A*x+G+F;
end