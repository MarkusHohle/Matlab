function [N_a, T] = poiss_stepper_fast(n_0, k)

N_a   = [1:n_0];
%N_a_s = sort(N_a, 'descend');


R   = rand(n_0,1);
Tau = - (1./(k*N_a)).*log(R');
T   = cumsum(Tau);

figure
stairs(T,N_a,'k-','LineWidth',3)
xlabel('time')
ylabel('location/state')