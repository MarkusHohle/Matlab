function poiss_stepper_manyTimes(Nrep,n_0,k)

N_a   = [1:n_0];
N_a_s = sort(N_a, 'descend');
LN    = length(N_a);

Tall  = zeros(Nrep,LN);

figure

for i=1:Nrep
    
    R   = rand(n_0,1);
    Tau = - (1./(k*N_a_s)).*log(R');
    T   = cumsum(Tau);
    
    Tall(i,:) = T; 
    
    stairs(T,N_a_s,'-','Color',[0.9 0.9 0.9])
    hold on
    
end

Tave = mean(Tall);
%Nave = mean(Nall');

stairs(Tave,N_a_s,'-','Color',[0 0 0],'LineWidth',3)
xlabel('time')
ylabel('location/state')