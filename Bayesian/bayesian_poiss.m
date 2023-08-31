function [q, Pu] = bayesian_poiss(data,options)
%function estimates the probability of a poissonian test and gives 1sigma 
%errorbars. 
%function uses two different priors: uniform and normal

%input: data the sequence of number of events of a poissonian process, 
% CI is the optional confidence interval, i. e. 0.68 is 1 sigma etc  


%example:
%generate random binomial numbers
%lambda = 2
%data   = poissrnd(lambda, 1, 10);
%do the BPE
% [x1, y1] = bayesian_poiss(data);
%
%data   = poissrnd(lambda, 1, 2);
%[x2, y2] = bayesian_poiss(data);
%
%with prior
%[x12, y12] = bayesian_poiss(data,Prior = [x1; y1]);
%
%other CI
%[x12, y12] = bayesian_poiss(data,Prior = [x1; y1], CI = 0.9);

arguments
       data double {mustBeNumeric}
       options.CI (1,1) {mustBeNumeric} = 0.68
       options.Prior {mustBeNumeric} = 1 %uniform prior
    end


CI = options.CI;

t    = length(data);
n    = sum(data);
qmax = (mean(data) + 1) *5;
dq   = 1/(100*t);
q    = eps:dq:qmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%uniform prior
%
%avoiding overflow for large N!
if n < 10
    Pu    = (((q*t).^n)./factorial(n)).*exp(-q*t);
else
    Pulog = n*log(q*t) - sum(log(1:n)) - q*t; %log P(lam|data) 
    Pu    = exp(Pulog); %P(lam|data) 
end

%multiplying pdf with prior
Prior    = options.Prior;
[nr, nc] = size(Prior);

if nr && nc > 1
    if nr <= nc 
      Prior = Prior';
    end
    x1              = Prior(:,1);
    y1              = Prior(:,2);
    yint            = interp1(x1,y1,q, 'spline');
else
    yint = Prior;
end

Pu              = Pu.*yint;
value           = trapz(q,Pu);
%normalisation
Pu              = Pu./value;
[maxPu, whereu] = max(Pu);

%determining errors (gauss approx of pdf)
pu             = q(whereu); %most likely value
sig_u          = sqrt(pu);

%determining 1sigma error by integral%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%positive direction
tot = 0;
i   = 0;
while tot < CI/2 && (i + whereu)<length(q)
    
    i   = i + 1;
    tot = trapz(q(whereu:whereu + i),Pu(whereu:whereu + i)); 
    
end

%negative direction
tot = 0;
j   = 0;
while tot < CI/2 && (whereu - j)> 1
    
    j   = j + 1;
    tot = trapz(q(whereu - j:whereu),Pu(whereu - j:whereu)); 
    
end

sig_pos = q(whereu + i);
sig_neg = q(whereu - j);

diff_pos = sig_pos - pu;
diff_neg = pu - sig_neg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rpu = round(pu*100)/100;
rdiff_pos = round(diff_pos*100)/100;
rdiff_neg = round(diff_neg*100)/100;

left  = q(whereu - j:whereu + i);

%right = q
Right = Pu(whereu - j:whereu + i);

figure
area(left, Right,'FaceColor',[0.8 0.8 0.8])
hold on
plot([q(whereu) q(whereu)],[0 max(Pu)*1.1],'r-')
hold on
plot(q,Pu,'k-')
box on
hold on
plot([pu+sig_u pu+sig_u],[0 max(maxPu)*1.1],'k--')
hold on
plot([pu-sig_u pu-sig_u],[0 max(maxPu)*1.1],'k--')
xlabel('q')
ylabel('P(q|{D},I)')
ylim([0 max(Pu)*1.1])
title(['q = ' num2str(rpu) '^{+' num2str(rdiff_pos) '}_{-' num2str(rdiff_neg) '}, ' num2str(CI*100) '% confidence' ])
legend([num2str(CI*100) '% confidence range'], 'most likely value','','1\sigma Gaussian approx (n = \infty)','')





