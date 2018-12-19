function f = High_Power10( gamma,p )
%High Power Approximation
f=exp(-(2-p)*gamma)/(2-p)-10*exp(-10*p*gamma)/p;
end