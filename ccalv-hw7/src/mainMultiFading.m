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
num_samples = 100;

fm_lower = 20; %Hz
fm_upper = 200; %Hz

%% 1. Generate Envelope Samples

% Each row is another sample
envelope_samples_lower = [];
envelope_samples_upper = [];

for s=1:num_samples
% Generate envelopes
    [time_lower, envelope_lower, waveform_time_lower] = rayleighFading(fm_lower, 8192);
    [time_upper, envelope_upper, waveform_time_upper] = rayleighFading(fm_upper, 8192);

    % Scale the envelopes by the RMS value
    envelope_lower = 10.*log10((10.^(envelope_lower./10))./rms(10.^(envelope_lower./10)));
    envelope_upper = 10.*log10((10.^(envelope_upper./10))./rms(10.^(envelope_upper./10)));

    % Append to total list
    envelope_samples_lower = [envelope_samples_lower; envelope_lower];
    envelope_samples_upper = [envelope_samples_upper; envelope_upper];
end

%% 2. Determine Values
clc
ref_lower = -10; %dB
ref_upper = 0; %dB

lcr_lower_ref_lower = 0;
lcr_lower_ref_upper = 0;
lcr_upper_ref_lower = 0;
lcr_upper_ref_upper = 0;

afd_lower_ref_lower = 0;
afd_lower_ref_upper = 0;
afd_upper_ref_lower = 0;
afd_upper_ref_upper = 0;

% Lower first
lower_crossings = 0;
upper_crossings = 0;
% Variables that answer the question, "Are we current fading?"
lower_fading_flag = 0;
lower_fading_times = [];
upper_fading_flag = 0;
upper_fading_times = [];
for s=1:num_samples
    last_level = envelope_samples_lower(s, 1);
    for p=2:(N+1)
        current_level = envelope_samples_lower(s, p);
        %-10dB Negative Crossing
        if (last_level > ref_lower) && (current_level < ref_lower)
            lower_fading_flag = 1;
            lower_fading_time = 0;
        %-10dB Positive Crossing
        elseif (last_level < ref_lower) && (current_level > ref_lower)
            lower_crossings = lower_crossings + 1;
            lower_fading_flag = 0;
            lower_fading_times = [lower_fading_times lower_fading_time];
        end
        if lower_fading_flag
            lower_fading_time = lower_fading_time + 1;
        end
        %0dB Negative Crossing
        if (last_level > ref_upper) && (current_level < ref_upper)
            upper_fading_flag = 1;
            upper_fading_time = 0;
        %0dB Positive Crossing
        elseif (last_level < ref_upper) && (current_level > ref_upper)
            upper_crossings = upper_crossings + 1;
            upper_fading_flag = 0;
            upper_fading_times = [upper_fading_times upper_fading_time];
        end
        if upper_fading_flag
            upper_fading_time = upper_fading_time + 1;
        end
        last_level = current_level;
    end
end
lcr_lower_ref_lower = (lower_crossings/num_samples)/waveform_time_lower;
lcr_lower_ref_upper = (upper_crossings/num_samples)/waveform_time_lower;
afd_lower_ref_lower = mean(lower_fading_times.*(waveform_time_lower/num_samples));
afd_lower_ref_upper = mean(upper_fading_times.*(waveform_time_lower/num_samples));
disp(strcat("Level Crossing Rate for ref_level=",num2str(ref_lower),"dB, doppler_max=",num2str(fm_lower), "Hz is ",num2str(lcr_lower_ref_lower)));
disp(strcat("Level Crossing Rate for ref_level=",num2str(ref_upper),"dB, doppler_max=",num2str(fm_lower), "Hz is ",num2str(lcr_lower_ref_upper)));
disp(strcat("Average Fade Duration for ref_level=",num2str(ref_lower),"dB, doppler_max=",num2str(fm_lower), "Hz is ",num2str(afd_lower_ref_lower)));
disp(strcat("Average Fade Duration for ref_level=",num2str(ref_upper),"dB, doppler_max=",num2str(fm_lower), "Hz is ",num2str(afd_lower_ref_upper)));


% Now upper
lower_crossings = 0;
upper_crossings = 0;
% Variables that answer the question, "Are we current fading?"
lower_fading_flag = 0;
lower_fading_times = [];
upper_fading_flag = 0;
upper_fading_times = [];
for s=1:num_samples
    last_level = envelope_samples_upper(s, 1);
    for p=2:(N+1)
        current_level = envelope_samples_upper(s, p);
        %-10dB Negative Crossing
        if (last_level > ref_lower) && (current_level < ref_lower)
            lower_fading_flag = 1;
            lower_fading_time = 0;
        %-10dB Positive Crossing
        elseif (last_level < ref_lower) && (current_level > ref_lower)
            lower_crossings = lower_crossings + 1;
            lower_fading_flag = 0;
            lower_fading_times = [lower_fading_times lower_fading_time];
        end
        if lower_fading_flag
            lower_fading_time = lower_fading_time + 1;
        end
        %0dB Negative Crossing
        if (last_level > ref_upper) && (current_level < ref_upper)
            upper_fading_flag = 1;
            upper_fading_time = 0;
        %0dB Positive Crossing
        elseif (last_level < ref_upper) && (current_level > ref_upper)
            upper_crossings = upper_crossings + 1;
            upper_fading_flag = 0;
            upper_fading_times = [upper_fading_times upper_fading_time];
        end
        if upper_fading_flag
            upper_fading_time = upper_fading_time + 1;
        end
        last_level = current_level;
    end
end
lcr_upper_ref_lower = (lower_crossings/num_samples)/waveform_time_upper;
lcr_upper_ref_upper = (upper_crossings/num_samples)/waveform_time_upper;
afd_upper_ref_lower = mean(lower_fading_times.*(waveform_time_upper/num_samples));
afd_upper_ref_upper = mean(upper_fading_times.*(waveform_time_upper/num_samples));
disp(strcat("Level Crossing Rate for ref_level=",num2str(ref_lower),"dB, doppler_max=",num2str(fm_upper), "Hz is ",num2str(lcr_upper_ref_lower)));
disp(strcat("Level Crossing Rate for ref_level=",num2str(ref_upper),"dB, doppler_max=",num2str(fm_upper), "Hz is ",num2str(lcr_upper_ref_upper)));
disp(strcat("Average Fade Duration for ref_level=",num2str(ref_lower),"dB, doppler_max=",num2str(fm_upper), "Hz is ",num2str(afd_upper_ref_lower)));
disp(strcat("Average Fade Duration for ref_level=",num2str(ref_upper),"dB, doppler_max=",num2str(fm_upper), "Hz is ",num2str(afd_upper_ref_upper)));


%% 3. BER Calculations
fm = fm_lower;
num_bits =  10000;
bit_sequence = [];
for bit=1:num_bits
    bit_sequence = [bit_sequence randi([0 1])];
end
    
symbol_sequence = []; % For 4-QAM
bit_pair = [];
second_bit_flag = 0;
for bit=1:num_bits
    bit_pair = [bit_pair bit_sequence(bit)];
    if second_bit_flag
        second_bit_flag = 0;
        if (bit_pair(1) == 0) && (bit_pair(2) == 0)
            symbol_sequence = [symbol_sequence 1*exp(j*(pi/4))];
        elseif (bit_pair(1) == 0) && (bit_pair(2) == 1)
            symbol_sequence = [symbol_sequence 1*exp(j*((3*pi)/4))];
        elseif (bit_pair(1) == 1) && (bit_pair(2) == 0)
            symbol_sequence = [symbol_sequence 1*exp(j*((5*pi)/4))];
        else
            symbol_sequence = [symbol_sequence 1*exp(j*((7*pi)/4))];
        end
        bit_pair = [];
    else
        second_bit_flag = 1;
    end
end

snrs = -20:10;
bers = [];
for snr=1:numel(snrs)
    h_of_n = [];
    [channel_time, channel_envelope, channel_waveform_time] = rayleighFading(fm, 8192);
    for envelope_point=1:numel(channel_envelope)
        h_of_n = [h_of_n (channel_envelope(envelope_point)*exp(-j*(rand()*2*pi)))];
    end
    channel_samples = [];
    for symbol=1:numel(symbol_sequence)
        channel_samples = [channel_samples (h_of_n(symbol)*symbol_sequence(symbol)...
            +(10^(normrnd(snrs(snr), 3)/10)+10^(normrnd(snrs(snr), 3)/10)*j))];
    end

    recovered_symbols = channel_samples./h_of_n(1:numel(channel_samples));
    recovered_bits = [];
    for symbol=1:numel(recovered_symbols)
        if (real(recovered_symbols(symbol)) > 0) && (imag(recovered_symbols(symbol)) >= 0)
            recovered_bits = [recovered_bits 0 0];
        elseif (real(recovered_symbols(symbol)) <= 0) && (imag(recovered_symbols(symbol)) > 0)
            recovered_bits = [recovered_bits 0 1];
        elseif (real(recovered_symbols(symbol)) < 0) && (imag(recovered_symbols(symbol)) <= 0)
            recovered_bits = [recovered_bits 1 0];
        elseif (real(recovered_symbols(symbol)) >= 0) && (imag(recovered_symbols(symbol)) < 0)
            recovered_bits = [recovered_bits 1 1];
        end
    end
    num_correct = sum(recovered_bits == bit_sequence);
    ber = (numel(bit_sequence) - num_correct)/numel(bit_sequence);
    bers = [bers ber];
end

plot(bers, flip(snrs));
title("SNR vs BER");
set(gca, 'XScale', 'log');
xlabel("BER");
ylabel("SNR");
ylim([min(snrs), max(snrs)]);