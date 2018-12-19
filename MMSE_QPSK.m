function mmse = MMSE_QPSK( rou )
%MMSE of QPSK
fun=@(x) exp(-x.^2-rou/2)./cosh(2*sqrt(rou/2).*x);
mmse=(1/sqrt(pi))*integral(fun,sqrt(rou/2)-5/sqrt(2),sqrt(rou/2)+5/sqrt(2));
%fun=@(x) 2./(1+exp(4*sqrt(rou/2).*x)).*exp(-(x-sqrt(rou/2)).^2)/sqrt(pi);
%mmse=integral(fun,-inf,inf);
end