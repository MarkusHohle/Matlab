function [rho, p, Pmeans] = mean_var_analysis_example
%function performs Bayesian Model Selection (Chapt 4.3, D. S. Sivia, Data 
%Analysis) with a KO vs WT example

%output:
%rho    = P(data comes from ONE population regarding mean and variance|data)/
%         P(data comes from TWO different populations regarding mean and 
%         variance|data)
%
%p      = ordinary p-val from two-sample, two-tailed, unpooled t-test 
%      (as comparison)
%
%Pmeans = P(mean_WT>mean_KO|data)


load LinaData.mat %mouse data

%let's take for example thymus vs weight
%TvsW = LinaData.thym_Ncells;
TvsW = LinaData.thym_vs_weight;
%TvsW = LinaData.TCell_DN;
WT   = TvsW{2,1};
KO   = TvsW{2,2};
lWT  = length(WT);
lKO  = length(KO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot the data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[valWT, whereWT] = hist(WT,5);
[valKO, whereKO] = hist(KO,5);

figure
bar(whereWT,valWT,1,'FaceColor','none','Edgecolor','r','LineWidth',2)
hold on
bar(whereKO,valKO,1,'FaceColor','none','Edgecolor','b','LineWidth',2)
legend('WT','KO')
xlabel('thymus/weight')
ylim([0 1.1*max([valWT,valKO])])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%actual analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All         = [WT',KO'];
lAll        = length(All);
%plausible prior constrains about means and var, here just the data range
mu_max_All  = max(All);
mu_min_All  = min(All);
std_max_All = max(All) - min(All);
std_min_All = 0; 

mu_max_WT  = max(WT);
mu_min_WT  = min(WT);
std_max_WT = max(WT) - min(WT);
std_min_WT = 0; 

mu_max_KO  = max(KO);
mu_min_KO  = min(KO);
std_max_KO = max(KO) - min(KO);
std_min_KO = 0; 

Nstep = 100;

dmu_all  = [mu_max_All  - mu_min_All]/Nstep;
dstd_all = [std_max_All - std_min_All]/Nstep;

dmu_WT  = [mu_max_WT  - mu_min_WT]/Nstep;
dstd_WT = [std_max_WT - std_min_WT]/Nstep;

dmu_KO  = [mu_max_KO  - mu_min_KO]/Nstep;
dstd_KO = [std_max_KO - std_min_KO]/Nstep;


%filling Z with loop for integration is way faster than meshgrid and direct
%matrix operations
Z = zeros(Nstep+1,Nstep+1);
x = [mu_min_All:dmu_all:mu_max_All]; 
y = [std_min_All:dstd_all:std_max_All];

for i=1:Nstep+1
    mu = x(i);
    for j=1:Nstep+1
          std = y(j);        
          
          Z(i,j) = exp(-0.5*sum(((All - mu)./(std + eps)).^2))./...
                   ((sqrt(2*pi*std^2)^lAll) + eps);                 
    end
end

%note: Integration is over sigma and mu and NOT over "All" --> pdf is not
%neccessarily normalized!!  --> thus, don't get puzzled when finding Iall >
%one
Iall = trapz(y,trapz(x,Z,2));


Z = zeros(Nstep+1,Nstep+1);
x = [mu_min_WT:dmu_WT:mu_max_WT]; 
y = [std_min_WT:dstd_WT:std_max_WT];

for i=1:Nstep+1
    mu = x(i);
    for j=1:Nstep+1
          std = y(j);        
          
          Z(i,j) = exp(-0.5*sum(((WT - mu)./(std + eps)).^2))./...
                   ((sqrt(2*pi*std^2)^lWT) + eps);                 
    end
end

IWT = trapz(y,trapz(x,Z,2));



Z = zeros(Nstep+1,Nstep+1);
x = [mu_min_KO:dmu_KO:mu_max_KO]; 
y = [std_min_KO:dstd_KO:std_max_KO];

for i=1:Nstep+1
    mu = x(i);
    for j=1:Nstep+1
          std = y(j);        
          
          Z(i,j) = exp(-0.5*sum(((KO - mu)./(std + eps)).^2))./...
                   ((sqrt(2*pi*std^2)^lKO) + eps);                 
    end
end

IKO = trapz(y,trapz(x,Z,2));


rho = (Iall/(IKO*IWT)) *...
      (1/((mu_max_All - mu_min_All)*(std_max_All - std_min_All)))*...
      (mu_max_WT - mu_min_WT)*(mu_max_KO - mu_min_KO)*...
      (std_max_WT - std_min_WT)*(std_max_KO - std_min_KO);
  
%to compare: t-test
[~,p] = ttest2(WT,KO,'Vartype','unequal');
    

%%%%%%%testing P(mu_1 > mu_2|Data)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mu_1 refers to WT

S_WT2  = var(WT);
S_KO2  = var(KO);

Stot2  = S_WT2/lWT + S_KO2/lKO;

z_hat  = mean(WT) - mean(KO);

prefac = 1/sqrt(2*pi*Stot2);

%zeta  = [0:inf];%difference between the means;
fun    = @(zeta, Stot2, z_hat) exp(-0.5 * ((zeta-z_hat).^2)./Stot2 );
Imeans = integral(@(zeta)fun(zeta, Stot2, z_hat),0,Inf);

Pmeans = prefac*Imeans;















