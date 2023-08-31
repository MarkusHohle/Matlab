function []=plotpowerspec(S,t)
%evaluating the signal with a fast fourier transformation (fft) and plots
%the power spectrum
%input: signal vector S
%       time vector t

%stepsize of the time vector
dt=abs(t(end)-t(end-1));

%sampling rate for the signal (large numer=good)
N=1e+5;

%matlab has an incluted fft function
fft_signal=fft(S,N);

%output from fft has to be modified (has something to do with complex 
%numbers. you can ask me, if you like to know more)
%note the denotation .^2 to square each element separatly. Also required
%while multiplying each element of the matrix A: A.^2 ~=A*A!!!!!
power=abs(fft_signal).^2 /N;

%to convert sampling steps into Hz (works only in this way if t is in sec)
F=[0:N/2-1]/(N*dt);

figure
semilogx(F(1:N/2),power(1:N/2),'k-')
xlabel('frequency [Hz]')
ylabel('power')
