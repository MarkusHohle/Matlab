function [D]=phage(t, G, rates)
%function provides the phage/bacteria/immune system ODEs to be solved with
%the ODE solver. See Leung & Weitz, 2016
%Input: G, a vector containing the actual concentrations of the
%       reactants/products P (phages), I (immune response) and B(bacteria)
%       t, the time vector
%Output: the actual lhs of the ODE
%note: units in ml, hour


B=G(1);
P=G(2);
I=G(3);


%reaction rates 
r       = rates(1); % growth rate of beacteria at low densities
beta    = rates(2); % burst size phages
phi     = rates(3); % adsorption rate phages
w       = rates(4); % decay rate phages
epsilon = rates(5); % killing rate param for immune response
alpha   = rates(6); % max growth rate of immune response

Kc      = rates(7); % carrying capacity of bacteria
KI      = rates(8); % carrying capacity immune response
KN      = rates(9); % bact conc, when immune response is half its max
KD      = rates(10);% bact conc, at which immune resp is half as effective 



%The model ODEs
dB = r*B*(1-B/Kc) - phi*B*P - epsilon*I*B/(1+B/KD);

dP = beta*phi*B*P - w*P;

dI = alpha*I*(1-I/KI)*B/(B+KN);

D = [dB; dP; dI];
