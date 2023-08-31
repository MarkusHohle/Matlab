function [D]=glycolysis(t, G)
%function provides the glycolysis ODEs to be solved with the ODE
%solver.
%Input: G, a vector containing the actual concentrations of the
%       reactants/products X (ADP) and Y(F6P)
%       t, the time vector
%Output: the actual lhs of the glycolysis ODE


X=G(1);
Y=G(2);

%reaction rates (limit cycle) 
%a=0.06;
%b=0.6;

%reaction rates (stationary state)  
a=0.05;
b=0.6;

%noise
a=a+0.5*rand;
b=b+0.5*rand;


%The glycolysis ODEs
dX = -X + a*Y + X*X*Y;

dY = b -a*Y - X*X*Y;

D = [dX; dY];
