function [B, P, I]=solve_phage(tin,START,rates)
%function solves the phage/bacteria/immune system ODEs using the matlab 
%ode45 solver (solves ODEs with the Dormand-Rince approach of fourth/fifth
%order Runge–Kutta method of adaptive step size)
%units in ml and hour


%%eg use%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%eg generating Fig 3c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t       = [0:1e-3:100];
%B0      = 1e3;
%P0      = 1e5;
%I0      = 1;

%start = [B0 P0 I0];

%rates = zeros(10,1);

%rates(1) = 1;       % growth rate of beacteria at low densities
%rates(2) = 100;     % burst size phages
%rates(3) = 5e-8;    % adsorption rate phages
%rates(4) = 1;       % decay rate phages
%rates(5) = 1e-6;    % killing rate param for immune response
%rates(6) = 1;       % max growth rate of immune response

%rates(7) = 1e9;     % carrying capacity of bacteria
%rates(8) = 1.5e6;   % carrying capacity immune response
%rates(9) = 1e4;     % bact conc, when immune response is half its max
%rates(10)= 1e8;     % bact conc, at which immune resp is half as effective 

%solve_phage(t, start, rates);

%%%%also interesting:
%rates(5) = 0; %no immune action on bacteria, Fig 3a
%rates(3) = 2e-10; %generates Figs 4 a/b (dashed) 
%solve_phage(t, start, rates);

%%%%%generates Figs 4 a/b (solid) 
%start(3) = rates(8); rates(5) = 1e-6;
%start(1) = (rates(7)-rates(10))/2 + sqrt((rates(7)-rates(10))^2/4 ...
%- rates(7)*rates(10)*rates(8)*rates(5)/rates(1));
%solve_phage(t, start, rates);

%%%%also interesting:
%rates(5) = 0; %no immune action on bacteria, Fig a
%rates(3) = 3e-11; %generates Figs 4 c/d (dashed) 
%solve_phage(t, start, rates);

%%%%%generates Figs 4 c/d (solid) 
%rates(5) = 1e-6;
%start(1) = (rates(7)-rates(10))/2 + sqrt((rates(7)-rates(10))^2/4 ...
%- rates(7)*rates(10)*rates(8)*rates(5)/rates(1));
%solve_phage(t, start, rates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tout, M]=ode45(@(tout, M)phage(tout, M, rates), tin, START, rates);

B = M(:,1);%bacteria conc
P = M(:,2);%phage conc
I = M(:,3);%immune resp

All = [B P I];


%plot temporal evolution of the concentrations of the reactants/products
figure
plot(tout,B,'Color',[0.1 0.8 0.1],'LineWidth',2)
hold on
plot(tout,P,'Color',[0.1 0.1 0.8],'LineWidth',2)
hold on
plot(tout,I,'k-','LineWidth',2)
xlabel('time [hrs]')
ylabel('concentration [ml^{-1}]')
legend('bacteria','phages','immune response')
set(gca,'Yscale','log')
ylim([0.7 1.7*max(All(:))])

%plot phase portrait and the last point as a red cross
figure
plot3(B,P,I,'k-','LineWidth',2)
hold on
plot3(B(end),P(end),I(end),'r+','MarkerSize',4)
xlabel('bacteria [ml^{-1}]')
ylabel('phages [ml^{-1}]')
zlabel('immune response [ml^{-1}]')
title('phase plot')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
set(gca,'Zscale','log')
xlim([0.1 1.1*max(B)])
ylim([0.1 1.1*max(P)])
zlim([0.1 1.1*max(I)])


box on




