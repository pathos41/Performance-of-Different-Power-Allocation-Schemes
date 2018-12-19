function rou = Bisection_64_QAM_23_new( a,b,tol,eta,gamma )
max=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
for k=1:1+max
    half=(a+b)/2;  %bisection
    %mmse1=MMSE_64_QAM_23_new(a);
    %mmse2=MMSE_64_QAM_23_new(b);
    mmse3=MMSE_64_QAM_23_new(half);
    
    %f1=mmse1-eta/gamma;
    %f2=mmse2-eta/gamma;
    f3=mmse3-eta/gamma;
    
    if f3==0   %Find the root
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