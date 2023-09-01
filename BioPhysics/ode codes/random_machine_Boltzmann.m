function random_machine_Boltzmann(N,M)

%function illustrates the increase of entropy by rolling ONE die among N 
%dice each of M time steps by starting with a homogenious set and evolving
%with a poissonian time keeper (Gillespie alg)


%N dice, all set to state "three" for t=1
Dice    = 3* ones(N,1);
T       = zeros(M,1);
Entropy = zeros(M,1);
Mean    = 3* ones(M,1);


%some plots:
[number , ~] = hist(Dice,6);
figure
subplot(2,2,1)
bar([0:5],number, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
title('histogram initial states','FontSize',12,'FontWeight','bold')
xlabel('micro state','FontSize',12,'FontWeight','bold')
ylabel('#','FontSize',12,'FontWeight','bold')
ylim([0 1.05*max(number)])
xlim([0.1 6.8])

for n=2:1:M
    
    %SETTING TIME
    %generate a random number between 0 and 1 (uniformly dist)
    r=rand; 
    tau=-log(r); 
    T(n)=T(n-1)+tau;
    
    %SETTING STATE of one randomly choosen die
    %choose die
    die = randi(N); %integer between 1 and N
    
    %choose new state
    state = randi(6); %integer between 1 and 6
    
    old_state = Dice(die);
    diff      = old_state - state;
        
    %reset old state
    Dice(die) = state;
    
    %here starts the conservation of energy
    test = 1;
    
    while test~=0
              
        %choose another die 
        die2 = randi(N);
        %get old state of the second die
        old_state2 = Dice(die2);
        
        new_state = old_state2 + diff;
        
        %check whether state makes sense      
        if new_state>0
            
            %reset old state of second die        
            Dice(die2) = new_state;
            test = 0;
        end
       
    end
    
    %calculating entropy of the system
    maxEn = max(Dice); %largest energy
    [number , ~] = hist(Dice,maxEn);
    p = number./N;
    entropy    =  -( p*(log(p+eps))' )* N;
    Entropy(n) = entropy;
    
    %calculating internal energy of the system
    m = mean(Dice);
    Mean(n) = m;

    
          
end

Emax    = N*(log(maxEn));

%some plots:
[number , ~] = hist(Dice, maxEn );
subplot(2,2,2)
bar([1:maxEn],number, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
title('histogram final states','FontSize',12,'FontWeight','bold')
xlabel('micro state','FontSize',12,'FontWeight','bold')
ylabel('#','FontSize',12,'FontWeight','bold')
ylim([0 1.05*max(number)])

subplot(2,2,3)
plot(T, Entropy,'r-x','LineWidth',3)
hold on
plot([T(1) T(end)],[Emax Emax],'k-')
xlabel('time','FontSize',12,'FontWeight','bold')
ylabel('total entropy','FontSize',12,'FontWeight','bold')
xlim([0 1.05*T(end)])

subplot(2,2,4)
plot(T, Mean,'r-x','LineWidth',3)
xlabel('time','FontSize',12,'FontWeight','bold')
ylabel({'internal energy (mean)', 'per particle'},'FontSize',12,'FontWeight','bold')
xlim([0 1.05*T(end)])














