% XXX.m
%     Cameron J. Calv
%     
% Demo.m
%     Kapil R. Dandekar
%     Xaime Rivas Rey
%	  ECE-T512 - Wireless Systems Matlab Courseware  
%% 0. Housekeeping and Set Important Values
clear all, close all, clc

iValue = 1;
jValue = 1;

letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G'];
inner_rad = 10000;
radius = (2*(inner_rad))/sqrt(3); %Meters
N = iValue^2 + iValue*jValue + jValue^2;
D = sqrt(3*N)*radius;
numberOfClusters = 1;
letterCountsPerFrame = [];

rxPowersB1 = [];
rxPowersB1Shad = [];
rxPowersB1ShadWorst95 = [];
powerCallMin = -88; %dBm

baseStationLabels = {};
num_rx_antennas = 4;
ant_space = 0.5;


frequency = 1.8*10^9;


% Determine bit and symbol sequence
fm = 20;
num_bits =  10000;
bit_sequence = [];
bit_sequence_2 = [];
for bit=1:num_bits
    bit_sequence = [bit_sequence randi([0 1])];
    bit_sequence_2 = [bit_sequence_2 randi([0 1])];
end

symbol_sequence = []; % For 4-QAM
symbol_sequence_2 = [];
bit_pair = [];
bit_pair_2 = [];
second_bit_flag = 0;
for bit=1:num_bits
    bit_pair = [bit_pair bit_sequence(bit)];
    bit_pair_2 = [bit_pair_2 bit_sequence(bit)];
    if second_bit_flag
        second_bit_flag = 0;
        % First symbol stream
        if (bit_pair(1) == 0) && (bit_pair(2) == 0)
            symbol_sequence = [symbol_sequence 1*exp(j*(pi/4))];
        elseif (bit_pair(1) == 0) && (bit_pair(2) == 1)
            symbol_sequence = [symbol_sequence 1*exp(j*((3*pi)/4))];
        elseif (bit_pair(1) == 1) && (bit_pair(2) == 0)
            symbol_sequence = [symbol_sequence 1*exp(j*((5*pi)/4))];
        else
            symbol_sequence = [symbol_sequence 1*exp(j*((7*pi)/4))];
        end
        % Second symbol stream
        if (bit_pair_2(1) == 0) && (bit_pair_2(2) == 0)
            symbol_sequence_2 = [symbol_sequence_2 1*exp(j*(pi/4))];
        elseif (bit_pair_2(1) == 0) && (bit_pair_2(2) == 1)
            symbol_sequence_2 = [symbol_sequence_2 1*exp(j*((3*pi)/4))];
        elseif (bit_pair_2(1) == 1) && (bit_pair_2(2) == 0)
            symbol_sequence_2 = [symbol_sequence_2 1*exp(j*((5*pi)/4))];
        else
            symbol_sequence_2 = [symbol_sequence_2 1*exp(j*((7*pi)/4))];
        end
        bit_pair = [];
        bit_pair_2 = [];
    else
        second_bit_flag = 1;
    end
end
    

%% 1.Create data structure for matlab movie
video_flag = 0; % 0 (NO) || 1 (YES)
distance = 2000; % Meters
speed = 50; % Miles per hour
metric_speed = speed / 2.237; % Meters per second
% If each frame represents a second
numFrames = 100;
if video_flag
	videoName = 'mimo_simulation';
	movieFrames = moviein(numFrames);
	frame_interf = [];
end

%% 2. Set location of the mobile user for each frame in the movie
%   (gives the illusion of user movement)
positionMax = 150;
startCenter = (inner_rad*(-0.5 + 1.7*j))/2;
endCenter = (inner_rad*(0.5 + 1.6*i))/2;
mobilePos = linspace(endCenter,0-startCenter, numFrames );
% mobilePos = ones(1, numFrames)*(1e2)*(3.7834 - 2*i);

snrs = [];
bers = [];
bers_sd = [];
bers_mrc = [];
bers_spatial = [];
%% 3. Create the animation
for index = 1:numFrames    % Draw each frame in the movie
    clf;
    hold on;
    axis off;
    
    % Draw the serving cell and label it
    drawCell( 0, radius, 'A' );
    clusterCenters = 0;
    axis equal;
    letterCounts = zeros(1, N);

    
    % Draw the mobile user(s) at the appropriate location
    plot( mobilePos(index), 'rx' );
    
     % Draw a line connecting the center (basestation) of the serving cell 
    %    and the mobile user
    [servingCellNumber, servingTierNumber, servingCellCenter] = findServingCell( mobilePos(index), clusterCenters);
    line( [real(mobilePos(index)) real(servingCellCenter)],...
        [imag(mobilePos(index)) imag(servingCellCenter)],...
            'Color','blue',...
            'LineStyle','-',...
            'LineWidth', 1);
    hold off;
    
    
    baseStationLabels(numel(baseStationLabels)+1, 1) = {convertStringsToChars(letters(servingCellNumber)+"_"+num2str(servingTierNumber))};

    
    % Determine receive power  
    rxPow = calcRXPower(mobilePos(index), clusterCenters, frequency, 'path_loss_exponent_2.9');
    rxPowersB1 = [rxPowersB1 rxPow]; 

    noise = mw2dbm(1.381*10^(-23)*290*1*1); % rxpow/Thermal Noise
    snr = rxPow - noise; 
    h_of_n = [];
    H_of_n = []; % Each row is each antenna channel realization
    H_of_n_spatial = []; % Each element is a 2x2 cell array, element corresponds to channel point
    [channel_time, channel_envelope, channel_waveform_time] = rayleighFading(fm, 8192);
    % SISO
    for envelope_point=1:numel(channel_envelope)
        h_of_n = [h_of_n (channel_envelope(envelope_point)*exp(-j*(rand()*2*pi)))];
    end
    % MIMO
    for ant=1:num_rx_antennas
        h_of_n_mimo = [];
        [channel_time, channel_envelope, channel_waveform_time] = rayleighFading(fm, 8192);
        for envelope_point=1:numel(channel_envelope)
            h_of_n_mimo = [h_of_n_mimo (channel_envelope(envelope_point)*exp(-j*(rand()*2*pi)))];
        end
        H_of_n = [H_of_n; h_of_n_mimo];
    end
    % MIMO Spatial Mux
    [channel_time_spat_1, channel_envelope_spat_1, channel_waveform_time_spat_1] = rayleighFading(fm, 8192);
    [channel_time_spat_2, channel_envelope_spat_2, channel_waveform_time_spat_2] = rayleighFading(fm, 8192);
    [channel_time_spat_3, channel_envelope_spat_3, channel_waveform_time_spat_3] = rayleighFading(fm, 8192);
    [channel_time_spat_4, channel_envelope_spat_4, channel_waveform_time_spat_4] = rayleighFading(fm, 8192);
    for envelope_point=1:8193
        block = [(channel_envelope_spat_1(envelope_point)*exp(-j*(rand()*2*pi))),...
            (channel_envelope_spat_2(envelope_point)*exp(-j*(rand()*2*pi)));...
            (channel_envelope_spat_3(envelope_point)*exp(-j*(rand()*2*pi))),...
            (channel_envelope_spat_4(envelope_point)*exp(-j*(rand()*2*pi)))];
        H_of_n_spatial = [H_of_n_spatial {block}];
    end
            
    channel_samples = []; 
    channel_samples_mimo = []; % Now is a column vector
    channel_samples_spatial = [];
    snr_of_sample = [];
    for symbol=1:numel(symbol_sequence)
        channel_samples = [channel_samples (dbm2mw(rxPow)*h_of_n(symbol)*symbol_sequence(symbol)...
            +(10^(normrnd(-snr, 3)/10)+10^(normrnd(-snr, 3)/10)*j))];
    end
    for symbol=1:numel(symbol_sequence)
        N_of_n = [];
        noise_of_n_dB = [];
        for ant=1:num_rx_antennas
            noise_real = normrnd(-snr, 3);
            noise_imag = normrnd(-snr, 3);
            noise_of_n_dB = [noise_of_n_dB; noise_real+noise_imag*j];
            N_of_n = [N_of_n; (10^(noise_real/10)+10^(noise_imag/10)*j)];
        end
        snr_of_sample = [snr_of_sample rxPow-noise_of_n_dB];
        channel_samples_mimo = [channel_samples_mimo (dbm2mw(rxPow)*(H_of_n(:, symbol).*eye(num_rx_antennas))*...
            ones(num_rx_antennas, 1)*symbol_sequence(symbol)...
            +(N_of_n))];
    end
    for symbol=1:numel(symbol_sequence)
        channel_samples_spatial = [channel_samples_spatial (dbm2mw(rxPow)*(H_of_n_spatial{symbol})*...
            [symbol_sequence(symbol); symbol_sequence_2(symbol)]...
            +([(10^(normrnd(-snr, 3)/10)+10^(normrnd(-snr, 3)/10)*j);...
            (10^(normrnd(-snr, 3)/10)+10^(normrnd(-snr, 3)/10)*j)]))];
    end
    
    recovered_symbols_vec = [];
    for symbol=1:numel(symbol_sequence)
        recovered_symbols_vec = [recovered_symbols_vec ...
            inv(H_of_n(:, symbol).*eye(num_rx_antennas))*channel_samples_mimo(:,symbol)];
    end
    % Maximal ratio combining symbol recovery (least squares technique)
    recovered_symbols_MRC = [];
    for symbol=1:numel(symbol_sequence)
        recovered_symbols_MRC = [recovered_symbols_MRC ...
            (H_of_n(:, symbol)'*...
            channel_samples_mimo(:, symbol))/...
            (H_of_n(:, symbol)'*H_of_n(:, symbol))];
    end
    recovered_symbols_spatial = [];
    for symbol=1:numel(symbol_sequence)
        recovered_symbols_spatial = [recovered_symbols_spatial inv(H_of_n_spatial{symbol})*...
            channel_samples_spatial(:, symbol)];
    end
    recovered_symbols = channel_samples./h_of_n(1:numel(channel_samples));
    recovered_bits = [];
    recovered_bits_highest_SNR = [];
    recovered_bits_max_ratio = [];
    recovered_bits_spatial = [];
    for symbol=1:numel(symbol_sequence)
        % Selection diversity symbol recovery
        [~, lowest_snr_idx] = max(snr_of_sample(:, symbol));
        chosen_symbol_highest_SNR = recovered_symbols_vec(lowest_snr_idx, symbol);
        
        % Maximal ratio combining symbol recovery
        average_noise = mean(10.^((rxPow - snr_of_sample(:, symbol))./10));
        chosen_symbol_max_ratio = recovered_symbols_MRC(symbol);
        
        % Recover bitstream for all cases
        if (real(recovered_symbols(symbol)) > 0) && (imag(recovered_symbols(symbol)) >= 0)
            recovered_bits = [recovered_bits 0 0];
        elseif (real(recovered_symbols(symbol)) <= 0) && (imag(recovered_symbols(symbol)) > 0)
            recovered_bits = [recovered_bits 0 1];
        elseif (real(recovered_symbols(symbol)) < 0) && (imag(recovered_symbols(symbol)) <= 0)
            recovered_bits = [recovered_bits 1 0];
        elseif (real(recovered_symbols(symbol)) >= 0) && (imag(recovered_symbols(symbol)) < 0)
            recovered_bits = [recovered_bits 1 1];
        end
        if (real(chosen_symbol_highest_SNR) > 0) && (imag(chosen_symbol_highest_SNR) >= 0)
            recovered_bits_highest_SNR = [recovered_bits_highest_SNR 0 0];
        elseif (real(chosen_symbol_highest_SNR) <= 0) && (imag(chosen_symbol_highest_SNR) > 0)
            recovered_bits_highest_SNR = [recovered_bits_highest_SNR 0 1];
        elseif (real(chosen_symbol_highest_SNR) < 0) && (imag(chosen_symbol_highest_SNR) <= 0)
            recovered_bits_highest_SNR = [recovered_bits_highest_SNR 1 0];
        elseif (real(chosen_symbol_highest_SNR) >= 0) && (imag(chosen_symbol_highest_SNR) < 0)
            recovered_bits_highest_SNR = [recovered_bits_highest_SNR 1 1];
        end
        if (real(chosen_symbol_max_ratio) > 0) && (imag(chosen_symbol_max_ratio) >= 0)
            recovered_bits_max_ratio = [recovered_bits_max_ratio 0 0];
        elseif (real(chosen_symbol_max_ratio) <= 0) && (imag(chosen_symbol_max_ratio) > 0)
            recovered_bits_max_ratio = [recovered_bits_max_ratio 0 1];
        elseif (real(chosen_symbol_max_ratio) < 0) && (imag(chosen_symbol_max_ratio) <= 0)
            recovered_bits_max_ratio = [recovered_bits_max_ratio 1 0];
        elseif (real(chosen_symbol_max_ratio) >= 0) && (imag(chosen_symbol_max_ratio) < 0)
            recovered_bits_max_ratio = [recovered_bits_max_ratio 1 1];
        end
        top_stream = [];
        bottom_stream = [];
        if (real(recovered_symbols_spatial(1, symbol)) > 0) && (imag(recovered_symbols_spatial(1, symbol)) >= 0)
            top_stream = [0 0];
        elseif (real(recovered_symbols_spatial(1, symbol)) <= 0) && (imag(recovered_symbols_spatial(1, symbol)) > 0)
            top_stream = [0 1];
        elseif (real(recovered_symbols_spatial(1, symbol)) < 0) && (imag(recovered_symbols_spatial(1, symbol)) <= 0)
            top_stream = [1 0];
        elseif (real(recovered_symbols_spatial(1, symbol)) >= 0) && (imag(recovered_symbols_spatial(1, symbol)) < 0)
            top_stream = [1 1];
        end
        if (real(recovered_symbols_spatial(2, symbol)) > 0) && (imag(recovered_symbols_spatial(2, symbol)) >= 0)
            bottom_stream = [0 0];
        elseif (real(recovered_symbols_spatial(2, symbol)) <= 0) && (imag(recovered_symbols_spatial(2, symbol)) > 0)
            bottom_stream = [0 1];
        elseif (real(recovered_symbols_spatial(2, symbol)) < 0) && (imag(recovered_symbols_spatial(2, symbol)) <= 0)
            bottom_stream = [1 0];
        elseif (real(recovered_symbols_spatial(2, symbol)) >= 0) && (imag(recovered_symbols_spatial(2, symbol)) < 0)
            bottom_stream = [1 1];
        end
        recovered_bits_spatial = [recovered_bits_spatial [top_stream; bottom_stream]];
    end
    num_correct = sum(recovered_bits == bit_sequence);
    num_correct_sd = sum(recovered_bits_highest_SNR == bit_sequence);
    num_correct_mrc = sum(recovered_bits_max_ratio == bit_sequence);
    num_correct_spatial_1 = sum(recovered_bits_spatial(1, :) == bit_sequence);
    num_correct_spatial_2 = sum(recovered_bits_spatial(2, :) == bit_sequence_2);
    ber = (numel(bit_sequence) - num_correct)/numel(bit_sequence);
    ber_sd = (numel(bit_sequence) - num_correct_sd)/numel(bit_sequence);
    ber_mrc = (numel(bit_sequence) - num_correct_mrc)/numel(bit_sequence);
    ber_spatial = ((numel(bit_sequence)+numel(bit_sequence_2)) -...
        (num_correct_spatial_1 + num_correct_spatial_2))/(numel(bit_sequence)+numel(bit_sequence_2));
    
    snrs = [snrs snr];
    bers = [bers ber];
    bers_sd = [bers_sd ber_sd];
    bers_mrc = [bers_mrc ber_mrc];
    bers_spatial = [bers_spatial ber_spatial];
    
    % Capture the frame for the movie
	movieFrames(:, index) = getframe(gcf);		
end

%% 4. Write the movie to a file (videoName.avi)
if video_flag
	v = VideoWriter(videoName);
	% Modify as needed to make your animation faster/slower
	v.FrameRate = 10;
	open(v);
	writeVideo(v,movieFrames);
	close(v);
end

%% 5. Plots
figure(2)
[sorted_snrs, sorted_idxs] = sort(snrs);
plot(sorted_snrs, bers(sorted_idxs));
hold on
plot(sorted_snrs, bers_sd(sorted_idxs));
plot(sorted_snrs, bers_mrc(sorted_idxs));
plot(sorted_snrs, bers_spatial(sorted_idxs));
title(strcat("BER vs SNR for Noise=", num2str(noise)));
legend("SISO", "Diversity Combining", "MRC", "2x2 Spatial Mux");
xlabel("SNR");
ylabel("BER");

figure(3)
plot(snrs);
title(strcat("SNR vs Moble Position for Noise=", num2str(noise)));
xlabel("Frame Number");
ylabel("SNR");

figure(4)
plot(bers);
hold on;
plot(bers_sd);
plot(bers_mrc);
plot(bers_spatial);
title("BER vs Moble Position");
legend("SISO", "Diversity Combining", "MRC", "2x2 Spatial Mux");
xlabel("Frame Number");
ylabel("BER");

