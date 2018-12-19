function rou = Bisection_QPSK( a,b,tol,eta,gamma )
max=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
for k=1:max+1
    half=(a+b)/2;  %bisection
    %mmse1=MMSE_QPSK(a);
    %mmse2=MMSE_QPSK(b);
    mmse3=MMSE_QPSK(half);
    
    %f1=mmse1-eta/gamma;
    %f2=mmse2-eta/gamma;
    f3=mmse3-eta/gamma;
    
    if f3==0   %Find the optimal point
        rou=half;
        break
    elseif f3<0
        b=half;
    else
        a=half;
    end
    if b-a<tol
        rou=half;
        break
    end
end
end
