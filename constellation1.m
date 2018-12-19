function p_k = constellation1( p,g_k,M_k,lamda )
%Constellation WF
if lamda<g_k*(1-1/M_k)
    p_k=1/(2*g_k*p)*(sqrt((M_k-1)^2+(4*g_k/lamda)*(M_k-1))-(M_k+1));
else
    p_k=0;
end
end