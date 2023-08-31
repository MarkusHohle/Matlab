function [X,Y]=solve_glycolysis(tin,START)
%function solves the glycolysis ODEs using the matlab ode45 solver 
%(solves ODEs with the Dormand-Rince approach of fourth/fifth order 
%Runge–Kutta method of adaptive step size)

%eg use solve_glycolysis([0:1e-3:100],[1 0]);
%note: timestep should be not too large in order to get a proper Fourier 
%analysis 

[tout, M]=ode45(@glycolysis, tin, START);

X=M(:,1);%[ADP]
Y=M(:,2);%[F6P]


%plot temporal evolution of the concentrations of the reactants/products
figure
plot(tout,X,'r-')
hold on
plot(tout,Y,'b-')
xlabel('time')
ylabel('concentration')
legend('ADP','F6P')

%plot X vs Y and the last point as a red cross
figure
plot(X,Y,'k-','LineWidth',2)
hold on
plot(X(end),Y(end),'r+','MarkerSize',4)
xlabel('ADP')
ylabel('F6P')
title('phase plot')

%Fourier analysis of the temporal evolution of X and Y
plotpowerspec(X,tout);
plotpowerspec(Y,tout);