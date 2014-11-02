clear all; close all;
load('HW5P1.mat')

%% Problem 1
% a)

plot(tax, oxr)

%%
% b)

sampling_interval = tax(2) - tax(1) % year^-1
Fs = 1 / sampling_interval;
Nyquist_frequency = 0.5 / sampling_interval

%%
% c)

ts = oxr - mean(oxr);
N = length(ts);
tsdft = fft(ts);
tsdft = tsdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(tsdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;

plot(freq, psdx)

%%
% d)

% 4000y
% 100,000y
% 2.6my

%%
% e)

w = [0; 1; 2; 2; 2; 2; 1; 0] / 10;
convts = conv(ts, w, 'same')

figure
plot(tax, ts, 'b')
hold on
plot(tax, convts, 'r')
hold off

%% Problem 2
% a)
load('HW5P2.mat')
N = size(sr, 1);
F = 100; % Hz
dt = 1/100;
t = dt:dt:N*dt;

plot(t, sr)

%% 
% b)

fs = 1/F;
[b, a] = butter(5, fs);
sr_lf = filter(b, a, sr);

plot(t, sr_lf)
