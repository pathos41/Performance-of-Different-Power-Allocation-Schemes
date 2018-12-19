%---------------------------Exact mercury/wf-------------------------------
tol=1e-5;   %Tolerance

c=0;   %Original lower bound for internal bisection
e=0;   %Original lower bound for internal bisection

d=500;  %Original upper bound for internal bisection
f=500;  %Original upper bound for internal bisection

Gamma1_dB=-20:20;
Gamma2_dB=Gamma1_dB+10;
Gamma3_dB=Gamma2_dB+10;

%Gamma11_dB=Gamma1_dB-3;
%Gamma21_dB=Gamma11_dB+3;
%Gamma31_dB=Gamma21_dB+3;  %xx-fold gain

%Gamma11=10.^(Gamma11_dB/10);
%Gamma21=10.^(Gamma21_dB/10);
%Gamma31=10.^(Gamma31_dB/10);

Gamma1=10.^(Gamma1_dB/10);
Gamma2=10.^(Gamma2_dB/10);
Gamma3=10.^(Gamma3_dB/10);

p1=zeros(1,length(Gamma1));
p2=zeros(1,length(Gamma1));
p3=zeros(1,length(Gamma1));

C_1=zeros(1,length(Gamma1));
C_2=zeros(1,length(Gamma1));
C_3=zeros(1,length(Gamma1));
C_4=zeros(1,length(Gamma1));
for n=1:length(Gamma1)
    gamma1=Gamma1(n);  %16-QAM
    gamma2=Gamma2(n);
    gamma3=Gamma3(n);
    
    a=gamma1*MMSE_16_QAM_23_new(300*gamma1/111);  %Original lower bound for external bisection
    b=gamma3*MMSE_16_QAM_23_new(300*gamma1/111);  %Original upper bound for external bisection
    
    %max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    %d=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    %f=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    
    for k=1:20
        eta=(a+b)/2;  %bisection
        
        %rou_1a=Bisection_16_QAM_23_new(c,d,tol,a,gamma1);
        %rou_2a=Bisection_16_QAM_23_new(e,f,tol,a,gamma2);
        %rou_3a=Bisection_16_QAM_23_new(e,f,tol,a,gamma3);
        
        %rou_1b=Bisection_16_QAM_23_new(c,d,tol,b,gamma1);
        %rou_2b=Bisection_16_QAM_23_new(e,f,tol,b,gamma2);
        %rou_3b=Bisection_16_QAM_23_new(e,f,tol,b,gamma3);
        
        rou_1e=Bisection_16_QAM_23_new(c,d,tol,eta,gamma1);
        rou_2e=Bisection_16_QAM_23_new(e,f,tol,eta,gamma2);
        rou_3e=Bisection_16_QAM_23_new(e,f,tol,eta,gamma3);
        
        %f_a=(1/(3*gamma1))*rou_1a+(1/(3*gamma2))*rou_2a+(1/(3*gamma3))*rou_3a-1;
        %f_b=(1/(3*gamma1))*rou_1b+(1/(3*gamma2))*rou_2b+(1/(3*gamma3))*rou_3b-1;
        f_e=(1/(3*gamma1))*rou_1e+(1/(3*gamma2))*rou_2e+(1/(3*gamma3))*rou_3e-1;
        
        if f_e==0   %Find the root
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            p3(n)=rou_3e/gamma3;
            break
        elseif f_e<0
            b=eta;
        else
            a=eta;
        end
        if abs(f_e)<1e-3
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            p3(n)=rou_3e/gamma3;
            break
        end
    end
end

%----------------------Constallation Constrained WF------------------------
p=10.^(Gamma1_dB/10);  %16-QAM
Pt=3;  %Total power

g1=1;
g2=10;
g3=100;  %Channel gain

M=16;  %Constellation order
tol=1e-5;   %Tolerance

p_1=zeros(1,length(p));
p_2=zeros(1,length(p));
p_3=zeros(1,length(p));
%-----------------------------Constellation WF-----------------------------
for n=1:length(p)
    a=0;   %Original lower bound
    b=g3*(1-1/M);   %Original upper bound
    max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    for k=1:max1+1
        lamda=(a+b)/2;   %Bisection
        %p_1a=constellation1(p(n),g1,M,a);   %Power allocation at a
        %p_2a=constellation1(p(n),g2,M,a);
        %p_3a=constellation1(p(n),g3,M,a);
        
        %p_1b=constellation1(p(n),g1,M,b);   %Power allocation at b
        %p_2b=constellation1(p(n),g2,M,b);
        %p_3b=constellation1(p(n),g3,M,b);
        
        p_1l=constellation1(p(n),g1,M,lamda);   %Power allocation at lamda
        p_2l=constellation1(p(n),g2,M,lamda);
        p_3l=constellation1(p(n),g3,M,lamda);
        
        %fa=p_1a+p_2a+p_3a-Pt;
        %fb=p_1b+p_2b+p_3b-Pt;
        fl=p_1l+p_2l+p_3l-Pt;
        
        if fl==0   %Find the root
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            p_3(n)=p_3l;
            break
        elseif fl<0
            b=lamda;
        else
            a=lamda;
        end
        if b-a<tol
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            p_3(n)=p_3l;
            break
        end
    end
end

%-----------------------------------p1-------------------------------------
C_2_1=zeros(length(Gamma1),1);
C_2_11=zeros(length(Gamma1),1);
C_2_12=zeros(length(Gamma1),1);
C_2_13=zeros(length(Gamma1),1);
C_2_14=zeros(length(Gamma1),1);

N=10000;
z=normrnd(0,1,N,1);  % Gaussian random variable

for n=1:length(Gamma1)
    sigma1=sqrt(1/Gamma1(n));
    %------------------------Exact mercury/wf------------------------------
    d_4_1=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(p1(n)/5);  % Distance matrix
    z1=sigma1*z;
    
    f1=zeros(length(d_4_1),1);
    b1=zeros(1,length(d_4_1));
    x1=zeros(N,1);
    
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_1)  % Number of columns
            for i=1:length(d_4_1)  % Number of rows
                f1(i)=exp(-d_4_1(i,j)^2/(2*sigma1^2)-z1(m)*d_4_1(i,j)/sigma1^2);
            end
            b1(j)=log(sum(f1))/log(2);  % First sum over i
        end
        x1(m)=sum(b1);  % Second sum over j
    end
    C_2_1(n)=4-0.5*mean(x1);  % Capacity of 16-QAM
    %------------------------Constellation wf------------------------------
    d_4_12=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(p_1(n)/5);  % Distance matrix
    
    f12=zeros(length(d_4_12),1);
    b12=zeros(1,length(d_4_12));
    x12=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_12)  % Number of columns
            for i=1:length(d_4_12)  % Number of rows
                f12(i)=exp(-d_4_12(i,j)^2/(2*sigma1^2)-z1(m)*d_4_12(i,j)/sigma1^2);
            end
            b12(j)=log(sum(f12))/log(2);  % First sum over i
        end
        x12(m)=sum(b12);  % Second sum over j
    end
    C_2_12(n)=4-0.5*mean(x12);  % Capacity of 16-QAM
    %--------------------Uniform power allocation--------------------------
    d_4_14=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(1/5);  % Distance matrix
    
    f14=zeros(length(d_4_14),1);
    b14=zeros(1,length(d_4_14));
    x14=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_14)  % Number of columns
            for i=1:length(d_4_14)  % Number of rows
                f14(i)=exp(-d_4_14(i,j)^2/(2*sigma1^2)-z1(m)*d_4_14(i,j)/sigma1^2);
            end
            b14(j)=log(sum(f14))/log(2);  % First sum over i
        end
        x14(m)=sum(b14);  % Second sum over j
    end
    C_2_14(n)=4-0.5*mean(x14);  % Capacity of 16-QAM
end

%-----------------------------------p2-------------------------------------
C_2_2=zeros(length(Gamma1),1);
C_2_21=zeros(length(Gamma1),1);
C_2_22=zeros(length(Gamma1),1);
C_2_24=zeros(length(Gamma1),1);
for n=1:length(Gamma1)
    sigma2=sqrt(1/Gamma2(n));
    %------------------------Exact mercury/wf------------------------------
    d_4_1=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(p2(n)/5);  % Distance matrix
    
    f1=zeros(length(d_4_1),1);
    b1=zeros(1,length(d_4_1));
    x1=zeros(N,1);
    
    z2=sigma2*z;  % Gaussian random variable
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_1)  % Number of columns
            for i=1:length(d_4_1)  % Number of rows
                f1(i)=exp(-d_4_1(i,j)^2/(2*sigma2^2)-z2(m)*d_4_1(i,j)/sigma2^2);
            end
            b1(j)=log(sum(f1))/log(2);  % First sum over i
        end
        x1(m)=sum(b1);  % Second sum over j
    end
    C_2_2(n)=4-0.5*mean(x1);  % Capacity of 16-QAM
    %------------------------Constellation wf------------------------------
    d_4_12=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(p_2(n)/5);  % Distance matrix
    
    f12=zeros(length(d_4_12),1);
    b12=zeros(1,length(d_4_12));
    x12=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_12)  % Number of columns
            for i=1:length(d_4_12)  % Number of rows
                f12(i)=exp(-d_4_12(i,j)^2/(2*sigma2^2)-z2(m)*d_4_12(i,j)/sigma2^2);
            end
            b12(j)=log(sum(f12))/log(2);  % First sum over i
        end
        x12(m)=sum(b12);  % Second sum over j
    end
    C_2_22(n)=4-0.5*mean(x12);  % Capacity of 16-QAM
    %--------------------Uniform power allocation--------------------------
    d_4_14=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(1/5);  % Distance matrix
    
    f14=zeros(length(d_4_14),1);
    b14=zeros(1,length(d_4_14));
    x14=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_14)  % Number of columns
            for i=1:length(d_4_14)  % Number of rows
                f14(i)=exp(-d_4_14(i,j)^2/(2*sigma2^2)-z2(m)*d_4_14(i,j)/sigma2^2);
            end
            b14(j)=log(sum(f14))/log(2);  % First sum over i
        end
        x14(m)=sum(b14);  % Second sum over j
    end
    C_2_24(n)=4-0.5*mean(x14);  % Capacity of 16-QAM
end
%-----------------------------------p3-------------------------------------
C_2_3=zeros(length(Gamma1),1);
C_2_31=zeros(length(Gamma1),1);
C_2_32=zeros(length(Gamma1),1);
C_2_34=zeros(length(Gamma1),1);
for n=1:length(Gamma1)
    sigma3=sqrt(1/Gamma3(n));
    %------------------------Exact mercury/wf------------------------------
    d_4_1=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(p3(n)/5);  % Distance matrix
    z3=sigma3*z;
    
    f1=zeros(length(d_4_1),1);
    b1=zeros(1,length(d_4_1));
    x1=zeros(N,1);
    
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_1)  % Number of columns
            for i=1:length(d_4_1)  % Number of rows
                f1(i)=exp(-d_4_1(i,j)^2/(2*sigma3^2)-z3(m)*d_4_1(i,j)/sigma3^2);
            end
            b1(j)=log(sum(f1))/log(2);  % First sum over i
        end
        x1(m)=sum(b1);  % Second sum over j
    end
    C_2_3(n)=4-0.5*mean(x1);  % Capacity of 16-QAM
    %------------------------Constellation wf------------------------------
    d_4_12=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(p_3(n)/5);  % Distance matrix
    
    f12=zeros(length(d_4_12),1);
    b12=zeros(1,length(d_4_12));
    x12=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_12)  % Number of columns
            for i=1:length(d_4_12)  % Number of rows
                f12(i)=exp(-d_4_12(i,j)^2/(2*sigma3^2)-z3(m)*d_4_12(i,j)/sigma3^2);
            end
            b12(j)=log(sum(f12))/log(2);  % First sum over i
        end
        x12(m)=sum(b12);  % Second sum over j
    end
    C_2_32(n)=4-0.5*mean(x12);  % Capacity of 16-QAM
    %------------------------Stronger channel------------------------------
    d_4_13=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(3/5);  % Distance matrix
    
    f13=zeros(length(d_4_13),1);
    b13=zeros(1,length(d_4_13));
    x13=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_13)  % Number of columns
            for i=1:length(d_4_13)  % Number of rows
                f13(i)=exp(-d_4_13(i,j)^2/(2*sigma3^2)-z3(m)*d_4_13(i,j)/sigma3^2);
            end
            b13(j)=log(sum(f13))/log(2);  % First sum over i
        end
        x13(m)=sum(b13);  % Second sum over j
    end
    C_2_13(n)=4-0.5*mean(x13);  % Capacity of 16-QAM
    %--------------------Uniform power allocation--------------------------
    d_4_14=[0:3;
        -1:2;
        -2:1;
        -3:0]'*2*sqrt(1/5);  % Distance matrix
    
    f14=zeros(length(d_4_14),1);
    b14=zeros(1,length(d_4_14));
    x14=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_4_14)  % Number of columns
            for i=1:length(d_4_14)  % Number of rows
                f14(i)=exp(-d_4_14(i,j)^2/(2*sigma3^2)-z3(m)*d_4_14(i,j)/sigma3^2);
            end
            b14(j)=log(sum(f14))/log(2);  % First sum over i
        end
        x14(m)=sum(b14);  % Second sum over j
    end
    C_2_34(n)=4-0.5*mean(x14);  % Capacity of 16-QAM
end
for n=1:length(Gamma1)
    C_1(n)=C_2_1(n)+C_2_2(n)+C_2_3(n);
    C_3(n)=C_2_12(n)+C_2_22(n)+C_2_32(n);
    C_4(n)=C_2_14(n)+C_2_24(n)+C_2_34(n);
end
plot(Gamma1_dB,C_1)
hold on
grid on
plot(Gamma1_dB,C_3,'--')
plot(Gamma1_dB,C_2_13,'*')
plot(Gamma1_dB,C_4,'-.')

xlabel('P/dB')
ylabel('Capacity')
legend('Exact mercury/waterfilling','AOPA','Stronger channel','Uniform power allocation')
