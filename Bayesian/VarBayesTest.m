function VarBayesTest(N, mu, sigma)
%Test function for the proof of principle applying variational bayes (VB)
%1)     data X={x1, x2, ...xN} are drawn from gaussian of random (given)
%       mu and sigma
%2)     estimating mu and sigma + their pdfs using VB (EM like algorithm)


X = normrnd(mu, sigma, 1, N);

xbar  = mean(X);
x2bar = mean(X.^2);

accur = 100;

%assuming mu was drawn from a gaussian and sigma from a gamma dist 
%(why? --> MaxEnt!)

%fixing hyperparams to small positive values (= largest ignorance to get 
%a broad prior)
lamb0 = eps;
mu0   = eps;
a0    = eps;
b0    = eps;
lambN = rand;
bN    = rand;

aN  = a0 + (N+1)/2;
muN = (lamb0*mu0 + N*xbar)/(lamb0 +N);

Tracker = zeros(1000,2);
Accur   = zeros(1000,1);
i=1;

while accur > 1e-7 && i < 1e+3
    
    %the approximate pdfs are only functions of aN, bN (Gamma) and muN,
    %lambN (Gaussian)
    Tracker(i,:) = [lambN bN];
    Accur(i)     = accur;
    i=i+1;
    
    %resetting in a circular(!!!) manner-----------------------------------
   
    bN_new  = b0 + 0.5*( (lamb0 + N)*(1/lambN + muN^2) - ...
        2*(lamb0*mu0 + N*xbar )*muN + N*x2bar +lamb0*mu0^2 );   
    
    lambN_new = (lamb0 + N)*aN/bN_new;
        
    accur = mean((abs(bN - bN_new)/abs(bN)) + ...
            (abs(lambN - lambN_new)/abs(lambN)));% + ...
        
    %overwriting for next iteration
    lambN = lambN_new;
    bN    = bN_new;

end

%final values
Tracker(i,:) = [lambN bN];
Accur(i)     = accur;

Tracker = Tracker(1:i,:);
Accur   = Accur(1:i);


%close all
figure
subplot(3,1,1)
plot([1:i],Tracker(:,1),'-x','LineWidth',3,'Color',[0.8 0.4 0])
xlabel('iteration')
ylabel('\lambda_N')
title('convergence check of parameters')


subplot(3,1,2)
plot([1:i],Tracker(:,2),'-x','LineWidth',3,'Color',[0.8 0.4 0])
xlabel('iteration')
ylabel('b_N')

subplot(3,1,3)
plot([1:i],100*Accur,'-x','LineWidth',3,'Color',[0.8 0.4 0])
set(gca,'YScale','log')
xlabel('iteration')
ylabel('mean accuracy [%]')
ylim([min(100*Accur)*0.8 max(100*Accur)*1.2])


%pdf of the precission tau (tau = 1/var) 
%creating the vector for the variance
Vardat = var(X);  % var of the data
dv = Vardat/N/10;
V  = [eps:dv:Vardat*2];

%calculating the pdf of V
gammaV = (1/gamma(aN))* (bN^aN)* ((1./V).^(aN-1)).*exp(-bN./V)  ;

%normalizing gamma (very crude, just for plotting). Filtering out NaNs that
%are caused by numerical inaccurracies
nn = isnan(gammaV);
f  = find(nn~=1);

gammaV = gammaV(f);
V      = V(f);

%actual normalization
Norm = trapz(V,gammaV);
gammaVN = gammaV/Norm;

[~, wo] = max(gammaVN);
%estimated variance from posterior approximation
VarEst = V(wo);

%plotting
figure
subplot(1,2,1)
plot(V,gammaVN,'-','LineWidth',2,'Color',[0.8 0.4 0])
hold on
plot([Vardat Vardat],[0 max(gammaVN)*1.05],'k-','LineWidth',2)
box on
ylim([0 max(gammaVN)*1.05])
xlim([V(1) V(end)])
xlabel('\sigma^2')
ylabel('P(\sigma^2|data)')
legend({'posterior for \sigma^2', ['variance of data\newline of size '...
    num2str(N)]},'FontWeight','bold')
title(['true variance: ' ...
    num2str(round((sigma^2)*100)/100) ',\newline estimated variance: ' ...
    num2str(round(VarEst*100)/100)])


%calculating the pdf of mu
Mudat = xbar; %mean of the data
dm    = abs(Mudat)/N/10;
Mu    = [min(X):dm:max(X)];

normM = normpdf(Mu,muN,sqrt(1/lambN));

[~, wo] = max(normM);
%estimated mean from posterior approximation
MEst = Mu(wo);

subplot(1,2,2)

if N<15
    
    [wert, wo] = hist(X);
    
else
    
    [wert, wo] = hist(X, round(N/4));
        
end

bar(wo, wert, 1, 'FaceColor', 'none')

%for plotting in the same range
normM = normM * max(wert) * 0.9;

hold on
plot(Mu,normM,'-','LineWidth',2,'Color',[0.8 0.4 0])
hold on
plot([Mudat Mudat],[0 max(normM)*1.05],'k-','LineWidth',2)
box on
ylim([0 max(wert)*1.05])
xlim([Mu(1) Mu(end)])
xlabel('\mu')
%ylabel('P(\mu|data)')
legend({'data' ,'P(\mu|data)', ['mean of data\newline of size '...
    num2str(N)]},'FontWeight','bold')
title(['true mean: ' ...
    num2str(round((mu)*100)/100) ',\newline estimated mean: ' ...
    num2str(round(MEst*100)/100)])





