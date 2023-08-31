function [L, W, E, T]=gillespieIV(L0,W0,E0,k1,k2,k3,k4,N)
%solves the stochastic predator-prey model (Lodka-Volterra/glycolysis)
%try for example [L W E T]=gillespieIV(200,5,200,0.1,0.025,0.05,0.5,1e5);


%allocating time vector
T=zeros(N,1);
%allocating vector for W, L and E
W=zeros(N,1);
L=zeros(N,1);
E=zeros(N,1);
%start with t=0
t=0;
Lt=L0;
Wt=W0;
Et=E0;


for i=1:1:N
    
    r1=rand;
    r2=rand;
    
    %total number of individuals (actually conserved, hence could be written outside the loop too)
    Ntot=Lt+Wt+Et;
    
    a1=k1*Lt*(Ntot-Lt-Wt)/Ntot/(Ntot-1);
    a2=k2*Wt/Ntot;
    a3=k3*Lt*Wt/Ntot/(Ntot-1);
    a4=k4*Lt*Wt/Ntot/(Ntot-1);
    
    a0=a1+a2+a3+a4;
    
    tau=1/(a0)*log(1/r1);
    
    if r2<a1/a0
        Et=Et-1;
        Lt=Lt+1;
    end
    
    if a1/a0<r2 && r2<(a1+a2)/a0
        Et=Et+1;
        Wt=Wt-1;
    end
    
    if (a1+a2)/a0<r2 && r2<(a1+a2+a3)/a0
        Lt=Lt-1;
        Wt=Wt+1;
    end
    
    if (a1+a2+a3)/a0<r2 && r2<1
        Lt=Lt-1;
        Et=Et+1;
    end
    
    t=t+tau;
    
    L(i)=Lt;
    E(i)=Et;
    W(i)=Wt;
    T(i)=t;
    
end

figure
stairs(T,L,'k-')
hold on
stairs(T,W,'r-')
%hold on
%stairs(T(f),E(f),'b-')
xlabel('time')
ylabel('number of molecules/individuals')
%legend('lambs', 'wolfs','empty')
legend('lambs', 'wolfs')


figure
plot(W,L,'k-')
xlabel('wolfs')
ylabel('lambs')
title('phase plot')


%Fourier analysis of the temporal evolution of X and Y
plotpowerspec(W,T);
plotpowerspec(L,T);

















