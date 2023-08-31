function poiss_stepper(n, nu)
%function simulates the poissonian stepper using the Gillespi Algorithm
%input:     - number of steps/states n
%           - hopping rate nu

%setting up a time and stae vector
N = zeros(n+1,1);
T = zeros(n+1,1);

for i = 2:n+1
    
    r = rand;
    tau = - (1/nu)*log(r);
    
    T(i)   = T(i-1) + tau;
    N(i) = N(i-1) + 1;
    
end

%figure
stairs(T,N,'k-','LineWidth',3)
xlabel('time')
ylabel('location/state')


%Homework: Write the code in a more efficient way, i. e. without using
%loops!