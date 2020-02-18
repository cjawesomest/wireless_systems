% main.m
%     Cameron J. Calv
%     
% Demo.m
%     Kapil R. Dandekar
%     Xaime Rivas Rey
%	  ECE-T512 - Wireless Systems Matlab Courseware
%        
%     This program illustrates how the Matlab courseware can be used to
% draw a cell and animate the trajectory of a mobile user moving across
% that cell. The program also shows how to record a video.
%% 0. Housekeeping and Set Important Values
clear all, close all, clc

N = 8192;

fm1 = 20; %Hz
fm2 = 200; %Hz
%% 1. Lower Max Shift
fm = fm1;
delta_f = (2*fm)/(N-1);
waveform_time = 1/(delta_f);

noise_freqs = linspace(-fm, fm, N+1);

% In_Phase
noise_sources = [];
for i=1:(N/2)
    noise_sources = [noise_sources (normrnd(0, 1)+j*normrnd(0, 1))];
end
noise_sources = [noise_sources 1 flip(conj(noise_sources))];

doppler_spectrum = noise_sources.*sqrt(1.5./(pi.*fm.*sqrt(1-((noise_freqs./fm).^2))));
last_slope = (doppler_spectrum(N)-doppler_spectrum(N-1))/(noise_freqs(N)-noise_freqs(N-1));
doppler_spectrum(N+1) = doppler_spectrum(N)+last_slope;
doppler_spectrum(1) = conj(doppler_spectrum(N+1));

in_phase = ifft(ifftshift(doppler_spectrum));

% Quadrature
noise_sources = [];
for i=1:(N/2)
    noise_sources = [noise_sources (normrnd(0, 1)+j*normrnd(0, 1))];
end
noise_sources = [noise_sources 1 flip(conj(noise_sources))];

doppler_spectrum = noise_sources.*sqrt(1.5./(pi.*fm.*sqrt(1-((noise_freqs./fm).^2))));
last_slope = (doppler_spectrum(N)-doppler_spectrum(N-1))/(noise_freqs(N)-noise_freqs(N-1));
doppler_spectrum(N+1) = doppler_spectrum(N)+last_slope;
doppler_spectrum(1) = conj(doppler_spectrum(N+1));

quadrature = ifft(ifftshift(doppler_spectrum));

time_series = sqrt((in_phase.^2) + (quadrature.^2));
subplot(2, 2, 1);
plot(linspace(0, waveform_time, numel(time_series)), 10*log10(time_series));
xlabel("Time (sec)");
xlim([0 6]);
ylabel("Relative RX Power (dB)");
title("Rayleigh Envelope Zoomed (f_{m}=20 Hz)");

subplot(2, 2, 3);
plot(linspace(0, waveform_time, numel(time_series)), 10*log10(time_series));
xlabel("Time (sec)");
xlim([0 waveform_time]);
ylabel("Relative RX Power (dB)");
title("Rayleigh Envelope (f_{m}=20 Hz)");

%% 2. Higher Max Shift
fm = fm2;
delta_f = (2*fm)/(N-1);
waveform_time = 1/(delta_f);

noise_freqs = linspace(-fm, fm, N+1);
% In_Phase
noise_sources = [];
for i=1:(N/2)
    noise_sources = [noise_sources (normrnd(0, 1)+j*normrnd(0, 1))];
end
noise_sources = [noise_sources 1 flip(conj(noise_sources))];

doppler_spectrum = noise_sources.*sqrt(1.5./(pi.*fm.*sqrt(1-((noise_freqs./fm).^2))));
last_slope = (doppler_spectrum(N)-doppler_spectrum(N-1))/(noise_freqs(N)-noise_freqs(N-1));
doppler_spectrum(N+1) = doppler_spectrum(N)+last_slope;
doppler_spectrum(1) = conj(doppler_spectrum(N+1));

in_phase = ifft(ifftshift(doppler_spectrum));

% Quadrature
noise_sources = [];
for i=1:(N/2)
    noise_sources = [noise_sources (normrnd(0, 1)+j*normrnd(0, 1))];
end
noise_sources = [noise_sources 1 flip(conj(noise_sources))];

doppler_spectrum = noise_sources.*sqrt(1.5./(pi.*fm.*sqrt(1-((noise_freqs./fm).^2))));
last_slope = (doppler_spectrum(N)-doppler_spectrum(N-1))/(noise_freqs(N)-noise_freqs(N-1));
doppler_spectrum(N+1) = doppler_spectrum(N)+last_slope;
doppler_spectrum(1) = conj(doppler_spectrum(N+1));

quadrature = ifft(ifftshift(doppler_spectrum));

subplot(2, 2, 2);
time_series = sqrt((in_phase.^2) + (quadrature.^2));
plot(linspace(0, waveform_time, numel(time_series)), 10*log10(time_series));
xlabel("Time (sec)");
xlim([0 0.5]);
ylabel("Relative RX Power (dB)");
title("Rayleigh Envelope Zoomed (f_{m}=200 Hz)");

subplot(2, 2, 4);
time_series = sqrt((in_phase.^2) + (quadrature.^2));
plot(linspace(0, waveform_time, numel(time_series)), 10*log10(time_series));
xlabel("Time (sec)");
xlim([0 waveform_time]);
ylabel("Relative RX Power (dB)");
title("Rayleigh Envelope (f_{m}=200 Hz)");


frequency = 2.4*10^9;