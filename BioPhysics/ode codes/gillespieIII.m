function [X Y T]=gillespieIII(X0,Y0,a,b,N)
%solves the stochastic Glycolysis ODEs
%reactants/products X (ADP), Y(F6P) and source term S0
%try eg [X Y T]=gillespieIII(20,20,0.01,0.75,3e3);

%allocating time vector
T=zeros(N,1);
%allocating vectors
X=zeros(N,1);
Y=zeros(N,1);
%start with t=0
t=0;
Xt=X0;
Yt=Y0;



for i=1:1:N
    
    r1=rand;
    r2=rand;
    
    a1=b;
    a2=a*Yt;
    a3=Yt*Xt^2;
    a4=Xt;
    
    a0=a1+a2+a3+a4;
        
    tau=1/(a0)*log(1/r1);
    
    if r2<a1/a0
        Yt=Yt+1;
    end
    
    if r2>a1/a0 && r2<(a2+a1)/a0
        Yt=Yt-1;
        Xt=Xt+1;
    end
    
    if r2>(a2+a1)/a0 && r2<(a1+a2+a3)/a0
        Yt=Yt-1;
        Xt=Xt+1;
    end
    
    if (a1+a2+a3)/a0<r2 %&& r2<1
        Xt=Xt-1;
    end
    
    if Xt<0
        Xt=0;
    end
    
    if Yt<0
        Yt=0;
    end
    
    
    t=t+tau;
    
    X(i)=Xt;
    Y(i)=Yt;
    T(i)=t;
    
end


figure
stairs(T,X,'k-')
hold on
stairs(T,Y,'r-')
xlabel('time')
ylabel('number of molecules/individuals')
legend('ADP(X)', 'F6P(Y)')

figure
plot(X,Y,'k-')
xlabel('ADP(X)')
ylabel('F6P(Y)')
title('phase plot')


%Fourier analysis of the temporal evolution of X and Y
plotpowerspec(X,T);
plotpowerspec(Y,T);

















