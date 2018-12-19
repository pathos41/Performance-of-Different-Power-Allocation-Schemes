function f = High_Power100( gamma,p )
%High Power Approximation
f=exp(-(2-p)*gamma)/(2-p)-100*exp(-100*p*gamma)/p;
end