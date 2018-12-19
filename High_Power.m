function f = High_Power( gamma,p )
%High Power Approximation
f=exp(-(2-p)*gamma)/(2-p)-2*exp(-2*p*gamma)/p;
end