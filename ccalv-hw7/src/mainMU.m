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
num_antennas = 10;
ant_space = 0.5;


frequency = 1.8*10^9;


% Determine bit and symbol sequence
fm = 20;
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
    

%% 1.Create data structure for matlab movie
video_flag = 0; % 0 (NO) || 1 (YES)
distance = 2000; % Meters
speed = 50; % Miles per hour
metric_speed = speed / 2.237; % Meters per second
% If each frame represents a second
numFrames = 100;
if video_flag
	videoName = 'mobile_user_vide';
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
    [channel_time, channel_envelope, channel_waveform_time] = rayleighFading(fm, 8192);
    for envelope_point=1:numel(channel_envelope)
        h_of_n = [h_of_n (channel_envelope(envelope_point)*exp(-j*(rand()*2*pi)))];
    end
    channel_samples = [];
    for symbol=1:numel(symbol_sequence)
        channel_samples = [channel_samples (dbm2mw(rxPow)*h_of_n(symbol)*symbol_sequence(symbol)...
            +(10^(normrnd(-snr, 3)/10)+10^(normrnd(-snr, 3)/10)*j))];
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
    
    snrs = [snrs snr];
    bers = [bers ber];
    
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
plot(snrs);
title(strcat("SNR vs Moble Position for Noise=", num2str(noise)));
xlabel("Frame Number");
ylabel("SNR");

figure(3)
plot(bers);
title("BER vs Moble Position");
xlabel("Frame Number");
ylabel("BER");

