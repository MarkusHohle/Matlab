function [D]=glycolysis_2(t, G, AB)
%function provides the glycolysis ODEs to be solved with the ODE
%solver.
%Input: G, a vector containing the actual concentrations of the
%       reactants/products X (ADP) and Y(F6P)
%       t, the time vector
%Output: the actual lhs of the glycolysis ODE


X=G(1);
Y=G(2);

%reaction rates 
a = AB(1);
b = AB(2);

%noise
%a=a+0.5*rand;
%b=b+0.5*rand;


%The glycolysis ODEs
dX = -X + a*Y + X*X*Y;

dY = b -a*Y - X*X*Y;

D = [dX; dY];
