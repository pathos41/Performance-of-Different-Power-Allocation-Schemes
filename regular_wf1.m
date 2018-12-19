function p_k = regular_wf1( p,g_k,lamda )
%Regular WF
if lamda<g_k*p
    p_k=1/lamda-1/(g_k*p);
else
    p_k=0;
end
end

