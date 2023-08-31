function [Evidence, Redchi2_min] = model_selection_example
%function performs Bayesian Model Selection (Chapt 4.1/4.2, D. S. Sivia, 
%Data Analysis) with an example of a clinical study (small n, low number of
%data points)

%output:
%
%Redchi2_min   = reduced chi2 for different models/same data from standart
%                least squares fit
%
%evidence      = (unnormalized) P(model|Data), priors P(model) are assumed
%                to be uniform

%close all

T    = readtable('tumor animal Fig2.txt');
conc = T.Var2; %ETS-TSP conc (mg/m3)
PA   = T.Var5; %tumor per animal
ePA  = T.Var7; %error tumor per animal


%%%%%%%%%%%%%%%%%%%%%%%%plot the data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure
%errorbar(conc(1:end-1),PA(1:end-1),ePA(1:end-1),'ko','MarkerFaceColor','k')
%hold on
%errorbar(conc(end),PA(end),ePA(end),'ko','MarkerFaceColor','r')
%xlim([-0.1 2.2])
%ylim([0 3])
%xlabel('ETS-TSP concentration [mg/m^3]','FontWeight','bold')
%ylabel('average tumors/animal','FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Model  = {'a*x+b','a*x.^2+b*x+c','a*x.^3+b*x.^2+c*x+d','A*exp(x*t)+C'};

%minimum and maximum expected params before(!) accessing the data (prior 
%information)  - a toy example!!
MinMax = {[0, 10; 0, 10], [-10 10; -10 10; -10 10], ...
    [-10 10; -10 10; -10 10; -10 10], [0 10; 0 10; 0 5]};

lM          = length(Model);
x           = conc(1:end-1);
y           = PA(1:end-1);
ey          = ePA(1:end-1);
n           = length(x);
Redchi2_min = zeros(lM,1);
Evidence    = zeros(lM,1);
Nstep = 1e+2; %for numerical integration

figure

for i=1:lM
    
    %actual fit
    current_model = Model{i};
    f             = fit(x, y, current_model);
    
    %evaluating the model symbolically for each x but unfixed params%%%%%%%
    model_f        = f(x);
    p              = length(coeffvalues(f));
    Redchi2_min(i) = sum(((model_f - y)./ey).^2) / (n - p);

    names         = coeffnames(f);
    syms(names)
    fun           = matlabFunction(eval(current_model));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %calculating integral of exp(-chi2/2) within MinMax of the params%%%%%%
    M     = MinMax{i};
    Ints  = zeros(p,Nstep+1);
    delta = zeros(p,1);
    
    %defining the integration area and increment
    for j=1:p
        
        delta(j)  = M(j,2) - M(j,1);
        dm        = delta(j)/Nstep;
        Ints(j,:) = [M(j,1):dm:M(j,2)];
       
    end
    
    switch i
        %------------------------------------------------------------------
        case 1
            Z = zeros(Nstep+1,Nstep+1);
            for aa = 1:Nstep
                a_model = Ints(1,aa);
                
                for bb = 1:Nstep
                    
                    b_model  = Ints(2,bb);
                    pred     = fun(a_model,b_model);
                    chi2     = sum(((pred - y)./ey).^2);
                    expchi   = exp(-chi2/2);
                    Z(aa,bb) = expchi;    
                    
                end
            end
            
            Grid{1}     = Ints(1,:); Grid{2} = Ints(2,:);
            integral    = trapezoidal_rule_nd_integral(Grid, Z, 2);
            Evidence(i) = (1/prod(delta)) * integral;

        %------------------------------------------------------------------
        case 2 
            Z = zeros(Nstep+1,Nstep+1,Nstep+1);
            for aa = 1:Nstep
                a_model = Ints(1,aa);
                
                for bb = 1:Nstep
                      b_model  = Ints(2,bb);

                    for cc = 1:Nstep
                    
                        c_model  = Ints(3,cc);
                        pred     = fun(a_model,b_model,c_model);
                        chi2     = sum(((pred - y)./ey).^2);
                        expchi   = exp(-chi2/2);
                        Z(aa,bb,cc) = expchi;    
                    end
                end
            end
            
            Grid{1} = Ints(1,:); Grid{2} = Ints(2,:);
            Grid{3} = Ints(3,:);
            integral    = trapezoidal_rule_nd_integral(Grid, Z, 3);
            Evidence(i) = (1/prod(delta)) * integral;

        %------------------------------------------------------------------
        case 3 
            Z = zeros(Nstep+1,Nstep+1,Nstep+1,Nstep+1);
            for aa = 1:Nstep
                a_model = Ints(1,aa);
                
                for bb = 1:Nstep
                      b_model  = Ints(2,bb);

                    for cc = 1:Nstep
                         c_model  = Ints(3,cc);
                         
                         for dd = 1:Nstep
                              d_model  = Ints(4,dd);

                            pred     = fun(a_model,b_model,c_model,...
                                       d_model);
                            chi2     = sum(((pred - y)./ey).^2);
                            expchi   = exp(-chi2/2);
                            Z(aa,bb,cc,dd) = expchi;   
                         end
                    end
                end
            end
            
            Grid{1} = Ints(1,:); Grid{2} = Ints(2,:);
            Grid{3} = Ints(3,:); Grid{4} = Ints(4,:);
            integral    = trapezoidal_rule_nd_integral(Grid, Z, 4);
            Evidence(i) = (1/prod(delta)) * integral;

        %------------------------------------------------------------------
        case 4
            Z = zeros(Nstep+1,Nstep+1,Nstep+1);
            for aa = 1:Nstep
                a_model = Ints(1,aa);
                
                for bb = 1:Nstep
                      b_model  = Ints(2,bb);

                    for cc = 1:Nstep
                    
                        c_model  = Ints(3,cc);
                        pred     = fun(a_model,b_model,c_model);
                        chi2     = sum(((pred - y)./ey).^2);
                        expchi   = exp(-chi2/2);
                        Z(aa,bb,cc) = expchi;    
                    end
                end
            end
            
            Grid{1} = Ints(1,:); Grid{2} = Ints(2,:);
            Grid{3} = Ints(3,:);
            integral    = trapezoidal_rule_nd_integral(Grid, Z, 3);
            Evidence(i) = (1/prod(delta)) * integral;
    end
    
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%plot the data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,i)
    errorbar(x,y,ey,'ko','MarkerFaceColor','k')
    hold on
    errorbar(conc(end),PA(end),ePA(end),'ko','MarkerFaceColor','r')
    hold on
    plot(f,'predobs',0.68)   
    xlim([-0.1 2.2])
    ylim([0 3])
    xlabel('ETS-TSP concentration [mg/m^3]','FontWeight','bold')
    ylabel('average tumors/animal','FontWeight','bold')  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
end




figure
bar([1:4],Evidence*1e+5,'FaceColor',[0.95 0.95 0.95])
hold on
plot([1:4],Redchi2_min,'kp--','MarkerFaceColor','k')
legend('P(model|{D},I)x10^5','\chi^2_{red}(min)')
set(gca,'Xtick',[1:4],'XTickLabel',Model,'XTickLabelRotation',45)
xlim([0 5])






