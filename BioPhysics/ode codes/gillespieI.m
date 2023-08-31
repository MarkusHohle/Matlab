function [A, T]=gillespieI(A0,k)
%simulates the reaction A-->B with rate k

%we need to define the vector A with a certain length
%but what is the length? a good approximation would be
%the number of initial molecules A0

A=zeros(A0,1)-1;
A(1,1)=A0;
T=zeros(A0,1);

t=0;
At=A0;
i=1;
while At>=0
    
        i=i+1;
        r=rand;
        tau=1/(At*k)*log(1/r);
        t=t+tau;
        At=At-1;
        
        A(i)=At;
        T(i)=t;
    
end

f=find(A>=0);

T=T(f);
A=A(f);

Aana=A0*exp(-k*T);

figure
stairs(T,A,'k-')
hold on
plot(T, Aana,'r-')
xlabel('time')
ylabel('number of molecules A')


