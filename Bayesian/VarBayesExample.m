function [Mean, Var] = VarBayesExample(data)
%function applies variational bayes (VB) in order to estimate the pdf of
%the mean and variance of the input data (matlab vetcor),
%(EM like algorithm, see HMM, the Rabiner paper, priors can be factorized 
%--> mean field approx)
%
%note: since only mean and variance are known --> assuming mu was drawn 
%from a gaussian and sigma from a gamma dist (why? --> MaxEnt!
%mean, var and normalization are the constrains for the Lagrangian)
%
%



%close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%uncommend, if Linas data wanted%%%%%%%%%%%%%
%load LinaData.mat %mouse data

%%let's take for example thymus vs weight
%TvsW = LinaData.thym_Ncells;
%TvsW = LinaData.thym_vs_weight;
%TvsW = LinaData.TCell_DN;
%WT   = TvsW{2,1};
%KO   = TvsW{2,2};

%data  = TvsW{2,2};

N      = length(data);
xbar   = mean(data);
x2bar  = mean(data.^2);
Vardat = var(data);

%X = normrnd(mu, sigma, 1, N);

accur = 400;



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

dv = Vardat/N/1000;
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
subplot(3,1,1)
plot(V,gammaVN,'-','LineWidth',2,'Color',[0.8 0.4 0])
hold on
plot([Vardat Vardat],[0 max(gammaVN)*1.05],'k-','LineWidth',2)
box on
%ylim([0 max(gammaVN)*1.05])
%xlim([V(1) V(end)])
xlabel('\sigma^2')
ylabel('P(\sigma^2|data)')
legend({'posterior for \sigma^2', ['variance of data\newline of size '...
    num2str(N)]},'FontWeight','bold')
title(['estimated variance: ' num2str(round(VarEst*1000)/1000)])


%calculating the pdf of mu
Mudat = xbar; %mean of the data
dm    = abs(Mudat)/N/1000;
Mu    = [min(data):dm:max(data)];

normM = normpdf(Mu,muN,sqrt(1/lambN));


[~, wo] = max(normM);
%estimated mean from posterior approximation
MEst = Mu(wo);

subplot(3,1,2)

if N<15
    
    [wert, wo] = hist(data);
    
else
    
    [wert, wo] = hist(data, round(N/4));
        
end

bar(wo, wert, 1, 'FaceColor', 'none')
ylim([0 max(wert)*1.05])
legend({'data'},'FontWeight','bold')
ylabel('#')

%for plotting in the same range
%normM = normM * max(wert) * 0.9;

subplot(3,1,3)
plot(Mu,normM,'-','LineWidth',2,'Color',[0.8 0.4 0])
hold on
plot([Mudat Mudat],[0 max(normM)*1.05],'k-','LineWidth',2)
box on
xlim([Mu(1) Mu(end)])
xlabel('\mu')
ylabel('P(\mu|data)')
legend({'posterior for \mu',['mean of data\newline of size ' num2str(N)]},...
    'FontWeight','bold')
title(['estimated mean: ' num2str(round(MEst*100)/100)])


Mean = [Mu' normM'];
Var  = [V' gammaVN'];


