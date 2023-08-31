function jumpGraphRandom(I)

A = [0 1 1 1 0; 1 0 1 0 0; 1 1 0 1 1; 1 0 1 0 0; 0 0 1 0 0];

coord = [2, 2; 4, 1; 3.8, 3; 1, 3.1; 4, 4];

[nr, nc] = size(A); 
R        = rand(nr, nc);
T        = R.*A; %creating random transition values

%leading to oscillations between 5 and 3
T(3,5)   = 100;
T(5,3)   = 100;

T        = T./sum(T,2);

e = rand(nr,1);  %creating random start probabilities
e = e/sum(e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nicer plot with size of edges corresponds to e and weigts are
%corresponding to color of the edges

%T is not neccessarily symmetric
Gupper = graph(T,'upper');
Glower = graph(T,'lower');

Eupper = Gupper.Edges;
Elower = Glower.Edges;

%creating gray scale color map from weights
Weights = (table2array(Eupper(:,2)) + table2array(Elower(:,2)))/2;
Weights = -Weights + max(Weights);
Weights = Weights./max(Weights);

Cmap   = repmat(Weights,1,3);

figure
T = sparse(T);
wgPlot(T,coord,'vertexWeight',e,'vertexScale',200,'edgeWidth',2,...
       'edgeColorMap', Cmap, 'vertexMetadata',e,'vertexColorMap',copper);
title('my graph')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


P = zeros(length(e),I+1);
P(:,1) = e/sum(e); %start values

for i = 1:I

    e        = T*e;
    P(:,i+1) = e/sum(e);

end

figure

for i = 1:5
    plot(1:I+1,P(i,:),'LineWidth',3)
    hold on
end

xlabel('iterations')
ylabel('probability')
legend({'1','2','3','4','5'})
hold off

