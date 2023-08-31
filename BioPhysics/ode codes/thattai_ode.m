function  thattai_ode ()
%function from erwin frey
% solves the deterministic single gene model of thattai and van oudenaarden

% parameters taken from the publication; time measured in units of seconds
b=2;        % b=2 or b=20
tau_r=120;  % mRNA half-time
tau_p=3600; % protein half-time

k1=0.01;
k2=log(2)/tau_r;
k3=b*k2;
k4=log(2)/tau_p;

% initial values
f0 = [0; 1200]; % play with the initial conditions
% time span
tmax = 200*60 % firts factor gives the timespan in units of minutes
tspan = [0; tmax];

% ode45 is the standard ODE-solver based on Runge-Kutta
fprintf('If you use octave, you must install odepkg manually. Download\n\thttp://sourceforge.net/projects/octave/files/Octave%%20Forge%%20Packages/Individual%%20Package%%20Releases/odepkg-0.8.0.tar.gz/download\nput the file in your current working directory and execute the following commands in octave:\n$ pkg prefix . .\n$ pkg install -local odepkg-0.8.0.tar.gz\nThe second command will take some time. Now you can use ode45 and related ode solving functions.\n\n');
options=[]; % time course plotted by default
[T,Y]=ode45(@f,tspan,f0,options,k1,k2,k3,k4);

% plot the result with time in minutes
subplot(3,1,1)
semilogx(T/60,Y(:,1),'o','color','r')
xlabel('Time (min)')
ylabel('[mRNA]')
xlim([T(1)/60  T(end)/60])
line(T/60,(k1/k2) .*(1-exp(-k2.*T)),'color','k')
subplot(3,1,2)
plot(T/60,Y(:,2),'o')
xlabel('Time (min)')
ylabel('[Protein]')

subplot(3,1,3)
% phase plane plotting in second figure
options = odeset('OutputFcn',@odephas2);
ode45(@f,tspan,f0,options,k1,k2,k3,k4);
xlabel('[mRNA]')
ylabel('[Protein]')


% -------------------------------------------------------------------------

function dydt = f(t,y,k1,k2,k3,k4)
dydt = [            k1-k2*y(1) 
                    k3*y(1)-k4*y(2)         ]; 
