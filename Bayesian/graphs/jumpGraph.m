function jumpGraph(I)

A = [0 1 1 1 0; 1 0 1 0 0; 1 1 0 1 1; 1 0 1 0 0; 0 0 1 0 0];
e = [1;1;1;1;1];

G = graph(A);

figure
plot(G)
title('my graph')

P = zeros(length(e),I+1);
P(:,1) = e/sum(e); %start values

for i = 1:I

    e        = A*e;
    P(:,i+1) = e/sum(e);

end

figure
plot(P','LineWidth',3)
xlabel('iterations')
ylabel('fraction of paths of length i')
legend({'1','2','3','4','5'})

