clc; clear all; close all;

%% Generating the room impulse response
fs = 16000;
M = fs/2 + 1;
framesize = 2048;

[B,A] = cheby2(4,20,[0.1 0.7]);
Hd = dfilt.df2t([zeros(1,6) B],A);

H = filter(Hd,log(0.99*rand(1,M)+0.01).*sign(randn(1,M))...
    .*exp(-0.002*(1:M)));
H = H/norm(H)*4;    
plot(0:1/fs:0.5,H);
xlabel('Time [sec]'); ylabel('Amplitude'); title('Room Impulse Response');
set(gcf, 'Color', [1 1 1])


%% Near end Speech Signal
load nearspeech

n = 1:length(v);
t = n/fs;
figure
plot(t,v);
axis([0 33.5 -1 1]); xlabel('Time [sec]'); ylabel('Amplitude');
title('Near-End Speech Signal');
set(gcf, 'Color', [1 1 1])


%% Far end Speech signal
load farspeech

farspeechEcho = filter(H,1,x);
figure
plot(t,farspeechEcho);
axis([0 33.5 -1 1]); xlabel('Time [sec]'); ylabel('Amplitude');
title('Far-End Echoed Speech Signal');
set(gcf, 'Color', [1 1 1])


%% Microphone signal
micSignal = farspeechEcho + v + 0.001*randn(length(v),1);
figure
plot(t,micSignal);
axis([0 33.5 -1 1]); xlabel('Time [sec]'); ylabel('Amplitude');
title('Microphone Signal');
set(gcf, 'Color', [1 1 1])


%% Adaptive Algorithm
mu = 0.025;
W0 = zeros(1,2048);
del = 0.01;
lam = 0.98;
x = x(1:length(W0)*floor(length(x)/length(W0)));
micSignal = micSignal(1:length(W0)*floor(length(micSignal)/length(W0)));


hFDAF = adaptfilt.fdaf(2048,mu,1,del,lam);
[y,e] = filter(hFDAF,x,micSignal);
n = 1:length(e);
t = n/fs;


%% ERLE Calculation
Hd2 = dfilt.dffir(ones(1,1000)); %returns a discrete-time, direct-form 
                                 %finite impulse response (FIR) filter, Hd2, 
                                 %with numerator coefficients, ones(1:1000)

erle = filter(Hd2,(e-v(1:length(e))).^2)./ ...
      (filter(Hd2,farspeechEcho(1:length(e)).^2));
erledB = -10*log10(erle);


%% Plotting Results
figure
subplot(411); plot(t,v(n),'g');
axis([0 33.5 -1 1]); ylabel('Amplitude'); title('Near-End Speech Signal');

subplot(412); plot(t,micSignal(n),'b'); 
axis([0 33.5 -1 1]); ylabel('Amplitude'); title('Microphone Signal');

subplot(413); plot(t,e(n),'r'); 
axis([0 33.5 -1 1]); xlabel('Time [sec]'); ylabel('Amplitude');
title('Output of Acoustic Echo Canceller');
set(gcf, 'Color', [1 1 1])

subplot(414); plot(t,erledB);
axis([0 33.5 0 70]); xlabel('Time [sec]'); ylabel('ERLE(dB)');
title('Echo Return Loss Enhancement');