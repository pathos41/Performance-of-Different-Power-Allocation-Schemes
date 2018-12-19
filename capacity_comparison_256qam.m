%---------------------------Exact mercury/wf-----------------------------
c=0;   %Original lower bound for internal bisection
e=0;   %Original lower bound for internal bisection

d=2500;  %Original upper bound for internal bisection
f=2500;  %Original upper bound for internal bisection

Gamma1_dB=-20:30;
Gamma2_dB=Gamma1_dB+3;   %xx-fold gain

%Gamma11_dB=Gamma1_dB-3;
%Gamma21_dB=Gamma11_dB+20;

%Gamma11=10.^(Gamma11_dB/10);
%Gamma21=10.^(Gamma21_dB/10);

Gamma1=10.^(Gamma1_dB/10);
Gamma2=10.^(Gamma2_dB/10);

p1=zeros(1,length(Gamma1));
p2=zeros(1,length(Gamma1));

p11=zeros(1,length(Gamma1));
p21=zeros(1,length(Gamma1));

C_1=zeros(1,length(Gamma1));
C_2=zeros(1,length(Gamma1));
C_3=zeros(1,length(Gamma1));
C_4=zeros(1,length(Gamma1));
C_5=zeros(1,length(Gamma1));

for n=1:length(Gamma1)
    gamma1=Gamma1(n);  %256-QAM
    gamma2=Gamma2(n);
    
    a=gamma1*MMSE_256_QAM_23_new(4*gamma1/3);  %Original lower bound for external bisection
    b=gamma2*MMSE_256_QAM_23_new(4*gamma1/3);  %Original upper bound for external bisection
    
    %max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    %d=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    %f=Bisection_4_PAM(0,100,1e-5,a,gamma2);  %Original upper bound
    
    for k=1:20
        eta=(a+b)/2;  %bisection
        
        %rou_1a=Bisection_256_QAM_23_new(c,d,1e-5,a,gamma1);
        %rou_2a=Bisection_256_QAM_23_new(e,f,1e-5,a,gamma2);
        
        %rou_1b=Bisection_256_QAM_23_new(c,d,1e-5,b,gamma1);
        %rou_2b=Bisection_256_QAM_23_new(e,f,1e-5,b,gamma2);
        
        rou_1e=Bisection_256_QAM_23_new(c,d,1e-5,eta,gamma1);
        rou_2e=Bisection_256_QAM_23_new(e,f,1e-5,eta,gamma2);
        
        %f_a=(1/(2*gamma1))*rou_1a+(1/(2*gamma2))*rou_2a-1;
        %f_b=(1/(2*gamma1))*rou_1b+(1/(2*gamma2))*rou_2b-1;
        f_e=(1/(2*gamma1))*rou_1e+(1/(2*gamma2))*rou_2e-1;
        
        if f_e==0   %Find the root
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            break
        elseif f_e<0
            b=eta;
        else
            a=eta;
        end
        if abs(f_e)<1e-3
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            break
        end
    end
end

%----------------------Constallation Constrained WF------------------------
p=10.^(Gamma1_dB/10);  %256-QAM
Pt=2;  %Total power

g1=1;
g2=2;  %Channel gain

M=256;  %Constellation order
tol=1e-5;   %Tolerance

p_1=zeros(1,length(p));
p_2=zeros(1,length(p));

%-----------------------------Constellation WF-----------------------------
for n=1:length(p)
    a=0;   %Original lower bound
    b=g2*(1-1/M);   %Original upper bound
    max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    for k=1:max1+1
        lamda=(a+b)/2;   %Bisection
        %p_1a=constellation1(p(n),g1,M,a);   %Power allocation at a
        %p_2a=constellation1(p(n),g2,M,a);
        
        %p_1b=constellation1(p(n),g1,M,b);   %Power allocation at b
        %p_2b=constellation1(p(n),g2,M,b);
        
        p_1l=constellation1(p(n),g1,M,lamda);   %Power allocation at lamda
        p_2l=constellation1(p(n),g2,M,lamda);
        
        %fa=p_1a+p_2a-Pt;
        %fb=p_1b+p_2b-Pt;
        fl=p_1l+p_2l-Pt;
        
        if fl==0   %Find the root
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            break
        elseif fl<0
            b=lamda;
        else
            a=lamda;
        end
        if b-a<tol
            p_1(n)=p_1l;
            p_2(n)=p_2l;
            break
        end
    end
end

%--------------------------------Regular WF--------------------------------
tol=1e-5;   %Tolerance
for n=1:length(p)
    c=0;   %Original lower bound
    d=g2;   %Original upper bound
    max2=-1+ceil((log(d-c)-log(tol))/log(2));  %Number of iterations
    for k=1:max2+1
        lamda_new=(c+d)/2;   %Bisection
        %p_1c=regular_wf1(p(n),g1,c);   %Power allocation at c
        %p_2c=regular_wf1(p(n),g2,c);
        
        %p_1d=regular_wf1(p(n),g1,d);   %Power allocation at d
        %p_2d=regular_wf1(p(n),g2,d);
        
        p_1l_new=regular_wf1(p(n),g1,lamda_new);   %Power allocation at lamda
        p_2l_new=regular_wf1(p(n),g2,lamda_new);
        
        %fc=p_1c+p_2c-Pt;   %Total power is 2
        %fd=p_1d+p_2d-Pt;
        fl_new=p_1l_new+p_2l_new-Pt;
        
        if fl_new==0   %Find the root
            p11(n)=p_1l_new;
            p21(n)=p_2l_new;
            break
        elseif fl_new<0
            d=lamda_new;
        else
            c=lamda_new;
        end
        if d-c<tol
            p11(n)=p_1l_new;
            p21(n)=p_2l_new;
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
C_2_15=zeros(length(Gamma1),1);

N=10000;
z=normrnd(0,1,N,1);  % Gaussian random variable

for n=1:length(Gamma1)
    sigma1=sqrt(1/Gamma1(n));
    %------------------------Exact mercury/wf------------------------------
    d_16_1=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(p1(n)/85);  % Distance matrix
    z1=sigma1*z;
    
    f1=zeros(length(d_16_1),1);
    b1=zeros(1,length(d_16_1));
    x1=zeros(N,1);
    
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_1)  % Number of columns
            for i=1:length(d_16_1)  % Number of rows
                f1(i)=exp(-d_16_1(i,j)^2/(2*sigma1^2)-z1(m)*d_16_1(i,j)/sigma1^2);
            end
            b1(j)=log(sum(f1))/log(2);  % First sum over i
        end
        x1(m)=sum(b1);  % Second sum over j
    end
    C_2_1(n)=8-0.125*mean(x1);  % Capacity of 256-QAM
    %------------------------Constellation wf------------------------------
    d_16_12=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(p_1(n)/85);  % Distance matrix
    
    f12=zeros(length(d_16_12),1);
    b12=zeros(1,length(d_16_12));
    x12=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_12)  % Number of columns
            for i=1:length(d_16_12)  % Number of rows
                f12(i)=exp(-d_16_12(i,j)^2/(2*sigma1^2)-z1(m)*d_16_12(i,j)/sigma1^2);
            end
            b12(j)=log(sum(f12))/log(2);  % First sum over i
        end
        x12(m)=sum(b12);  % Second sum over j
    end
    C_2_12(n)=8-0.125*mean(x12);  % Capacity of 256-QAM
    %--------------------Uniform power allocation--------------------------
    d_16_14=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(1/85);  % Distance matrix
    
    f14=zeros(length(d_16_14),1);
    b14=zeros(1,length(d_16_14));
    x14=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_14)  % Number of columns
            for i=1:length(d_16_14)  % Number of rows
                f14(i)=exp(-d_16_14(i,j)^2/(2*sigma1^2)-z1(m)*d_16_14(i,j)/sigma1^2);
            end
            b14(j)=log(sum(f14))/log(2);  % First sum over i
        end
        x14(m)=sum(b14);  % Second sum over j
    end
    C_2_14(n)=8-0.125*mean(x14);  % Capacity of 256-QAM
    %-----------------------------Regular wf-----------------------------------
    d_16_15=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(p11(n)/85);  % Distance matrix
    
    f15=zeros(length(d_16_15),1);
    b15=zeros(1,length(d_16_15));
    x15=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_15)  % Number of columns
            for i=1:length(d_16_15)  % Number of rows
                f15(i)=exp(-d_16_15(i,j)^2/(2*sigma1^2)-z1(m)*d_16_15(i,j)/sigma1^2);
            end
            b15(j)=log(sum(f15))/log(2);  % First sum over i
        end
        x15(m)=sum(b15);  % Second sum over j
    end
    C_2_15(n)=8-0.125*mean(x15);  % Capacity of 256-QAM
end

%-----------------------------------p2-------------------------------------
C_2_2=zeros(length(Gamma1),1);
C_2_21=zeros(length(Gamma1),1);
C_2_22=zeros(length(Gamma1),1);
C_2_24=zeros(length(Gamma1),1);
C_2_25=zeros(length(Gamma1),1);
for n=1:length(Gamma1)
    sigma2=sqrt(1/Gamma2(n));
    %------------------------Exact mercury/wf------------------------------
    d_16_1=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(p2(n)/85);  % Distance matrix
    
    f1=zeros(length(d_16_1),1);
    b1=zeros(1,length(d_16_1));
    x1=zeros(N,1);
    
    z2=sigma2*z;  % Gaussian random variable
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_1)  % Number of columns
            for i=1:length(d_16_1)  % Number of rows
                f1(i)=exp(-d_16_1(i,j)^2/(2*sigma2^2)-z2(m)*d_16_1(i,j)/sigma2^2);
            end
            b1(j)=log(sum(f1))/log(2);  % First sum over i
        end
        x1(m)=sum(b1);  % Second sum over j
    end
    C_2_2(n)=8-0.125*mean(x1);  % Capacity of 256-QAM
    %------------------------Constellation wf------------------------------
    d_16_12=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(p_2(n)/85);  % Distance matrix
    
    f12=zeros(length(d_16_12),1);
    b12=zeros(1,length(d_16_12));
    x12=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_12)  % Number of columns
            for i=1:length(d_16_12)  % Number of rows
                f12(i)=exp(-d_16_12(i,j)^2/(2*sigma2^2)-z2(m)*d_16_12(i,j)/sigma2^2);
            end
            b12(j)=log(sum(f12))/log(2);  % First sum over i
        end
        x12(m)=sum(b12);  % Second sum over j
    end
    C_2_22(n)=8-0.125*mean(x12);  % Capacity of 256-QAM
    %------------------------Stronger channel------------------------------
    d_16_13=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(2/85);  % Distance matrix
    
    f13=zeros(length(d_16_13),1);
    b13=zeros(1,length(d_16_13));
    x13=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_13)  % Number of columns
            for i=1:length(d_16_13)  % Number of rows
                f13(i)=exp(-d_16_13(i,j)^2/(2*sigma2^2)-z2(m)*d_16_13(i,j)/sigma2^2);
            end
            b13(j)=log(sum(f13))/log(2);  % First sum over i
        end
        x13(m)=sum(b13);  % Second sum over j
    end
    C_2_13(n)=8-0.125*mean(x13);  % Capacity of 256-QAM
    %--------------------Uniform power allocation--------------------------
    d_16_14=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(1/85);  % Distance matrix
    
    f14=zeros(length(d_16_14),1);
    b14=zeros(1,length(d_16_14));
    x14=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_14)  % Number of columns
            for i=1:length(d_16_14)  % Number of rows
                f14(i)=exp(-d_16_14(i,j)^2/(2*sigma2^2)-z2(m)*d_16_14(i,j)/sigma2^2);
            end
            b14(j)=log(sum(f14))/log(2);  % First sum over i
        end
        x14(m)=sum(b14);  % Second sum over j
    end
    C_2_24(n)=8-0.125*mean(x14);  % Capacity of 256-QAM
    %-----------------------------Regular wf-----------------------------------
    d_16_15=[0:15;
        -1:14;
        -2:13;
        -3:12;
        -4:11;
        -5:10;
        -6:9;
        -7:8;
        -8:7;
        -9:6;
        -10:5;
        -11:4;
        -12:3;
        -13:2;
        -14:1;
        -15:0;]'*2*sqrt(p21(n)/85);  % Distance matrix
    
    f15=zeros(length(d_16_15),1);
    b15=zeros(1,length(d_16_15));
    x15=zeros(N,1);
    
    for m=1:N  % Number of pages
        for j=1:length(d_16_15)  % Number of columns
            for i=1:length(d_16_15)  % Number of rows
                f15(i)=exp(-d_16_15(i,j)^2/(2*sigma2^2)-z2(m)*d_16_15(i,j)/sigma2^2);
            end
            b15(j)=log(sum(f15))/log(2);  % First sum over i
        end
        x15(m)=sum(b15);  % Second sum over j
    end
    C_2_25(n)=8-0.125*mean(x15);  % Capacity of 256-QAM
end

for n=1:length(Gamma1)
    C_1(n)=C_2_1(n)+C_2_2(n);
    C_3(n)=C_2_12(n)+C_2_22(n);
    C_4(n)=C_2_14(n)+C_2_24(n);
    C_5(n)=C_2_15(n)+C_2_25(n);
end
plot(Gamma1_dB,C_1)
hold on
grid on
plot(Gamma1_dB,C_3,'--')
plot(Gamma1_dB,C_2_13,'*')
plot(Gamma1_dB,C_4,'-.')
plot(Gamma1_dB,C_5,'.')

xlabel('P/dB')
ylabel('Capacity')
legend('Exact mercury/waterfilling','AOPA','Stronger channel','Uniform power allocation','Regular WF')
