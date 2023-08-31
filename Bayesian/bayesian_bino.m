function [q, Pu] = bayesian_bino(n,k,options)
%function estimates the probability of a binomial test and gives 1sigma 
%errorbars. 
%function uses two different priors: uniform and normal

%input: n length of the sequence; k number of successes, CI is the optional
%confidence interval, i. e. 0.68 is 1 sigma etc, Prior is a vector 
%containing x and y values of previous normalized BPE  


%example:
%generate random binomial numbers
%p = 0.2;
%n = 20;
%k = binornd(n,p);
%do the BPE
%[x1, y1] = bayesian_bino(n,k);
%
%n = 5;
%[x2, y2] = bayesian_bino(n,k);
%
%with prior
%[x12, y12] = bayesian_bino(n,k,Prior = [x1; y1]);
%
%other CI
%[x12, y12] = bayesian_bino(n,k,Prior = [x1; y1], CI = 0.9);

arguments
       n double
       k double
       options.CI (1,1) {mustBeNumeric} = 0.68
       options.Prior {mustBeNumeric} = 1 %uniform prior
    end


CI = options.CI;


dq = 1/(100*n);
q  = 0:dq:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pu = (q.^k).*(1 - q).^(n - k);

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
sig_u          = sqrt(pu*(1-pu)/n);

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





