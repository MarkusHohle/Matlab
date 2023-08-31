function [A, B, T]=gillespieIIb(A0,B0,k1,k2,N)
%solves reaction 
%   k1
%A <--> B
%   k2
%stochasticaly
%Input: A0=A(t=0) and B0=B(t=0)
%       reaktion rates k1 and k2
%       number of steps N

%[A B T]=gillespieII(22,15,0.1,0.5,200); reproduces the plot in the paper,
%since they set B0=1 implicitly


%allocating time vector
T=zeros(N,1);
%allocating vector for A and B
A=zeros(N,1);
B=zeros(N,1);
%start with t=0
t=0;
At=A0;
Bt=B0;

for i=1:1:N
    
    r1=rand;
    r2=rand;
    
    a0=At*k1 + Bt*k2;
    
    tau=1/(a0)*log(1/r1);
    
    t=t+tau;
    
    if r2<(k2*Bt)/a0
        At=At+1;
        %in the paper (Erban et al., A Practical guide to stochastic 
        %simulations of reaction-diffusion processes) B seems to has an 
        %infinit reservoir, so that At reaches steady state. If you don't 
        %like this, just uncommend the next lines containing Bt
        
        Bt=Bt-1;
          
    else
        At=At-1;
        Bt=Bt+1;
    end
        
    A(i)=At;
    B(i)=Bt;
    T(i)=t;
    
end


figure
stairs(T,A,'k-')
hold on
plot(T,B,'r-')
xlabel('time')
ylabel('number of molecules')
legend('A', 'B')

%figure
%plot(A(f),B(f),'k-')
%xlabel('A')
%ylabel('B')
