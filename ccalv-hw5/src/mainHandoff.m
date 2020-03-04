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
jValue = 2;

letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G'];
radius = (2*(1000))/sqrt(3); %Meters
N = iValue^2 + iValue*jValue + jValue^2;
D = sqrt(3*N)*radius;
numberOfClusters = 1;
letterCountsPerFrame = [];

rxPowersB1 = [];
rxPowersB2 = [];
rxPowersB1Shad = [];
rxPowersB2Shad = [];
rxPowersB1ShadWorst95 = [];
rxPowersB2ShadWorst95 = [];
powerCallMin = -88; %dBm

baseStationLabels = {};
frequency = 2.4*10^9;
    

%% 1.Create data structure for matlab movie
video_flag = 0; % 0 (NO) || 1 (YES)
distance = 2000; % Meters
speed = 50; % Miles per hour
metric_speed = speed / 2.237; % Meters per second
% If each frame represents a second
numFrames = ceil(distance/metric_speed);
if video_flag
	videoName = 'handoff_video';
	movieFrames = moviein(numFrames);
	frame_interf = [];
end

% Note that the coordinate system is x+j*y
%  i.e. X is real axis, and Y is imaginary axis
%  Remember to be consistent with your units

%% 2. Set location of the mobile user for each frame in the movie
%   (gives the illusion of user movement)
positionMax = 150;
startCenter = 0 - 0*j;
endCenter = 1.0e+03*(1.7321 + 1.0000*j);
mobilePos = linspace( 0-startCenter, endCenter, numFrames );

%For Homework 3
%Let's simulate handoff!


%% 3. Create the animation
for index = 1:numFrames    % Draw each frame in the movie
    clf;
    hold on;
    axis off;
    
    % Draw the serving cell and label it
%     drawCell( 0, 100, 'A_1' );
    clusterCenters = drawManyClusters(0, iValue, jValue, radius, numberOfClusters);
    letterCounts = zeros(1, N);
    
    % Draw the mobile user(s) at the appropriate location
    plot( mobilePos(index), 'x' );
    
   % Draw a line showing trajectory
    line( [real(startCenter) real(mobilePos(index))],...
        [imag(startCenter) imag(mobilePos(index))],...
        'Color', 'red');
     % Draw a line connecting the center (basestation) of the serving cell 
    %    and the mobile user
    [servingCellNumber, servingTierNumber, servingCellCenter] = findServingCell( mobilePos(index), clusterCenters);
    line( [real(mobilePos(index)) real(endCenter)],...
        [imag(mobilePos(index)) imag(endCenter)],...
            'Color','blue',...
            'LineStyle','-',...
            'LineWidth', 1);
    hold off;
    
    
    baseStationLabels(numel(baseStationLabels)+1, 1) = {convertStringsToChars(letters(servingCellNumber)+"_"+num2str(servingTierNumber))};

    
    % Determine receive power  
    rxPow = calcRXPower(mobilePos(index), startCenter, frequency, 'path_loss_exponent_2.9');
    rxPowersB1 = [rxPowersB1 rxPow];
    rxPow = calcRXPower(mobilePos(index), endCenter, frequency, 'path_loss_exponent_2.9');
    rxPowersB2 = [rxPowersB2 rxPow];
    rxPow = calcRXPower(mobilePos(index), startCenter, frequency, 'path_loss_exponent_2.9_shadowing');
    rxPowersB1Shad = [rxPowersB1Shad rxPow];
    rxPow = calcRXPower(mobilePos(index), endCenter, frequency, 'path_loss_exponent_2.9_shadowing');
    rxPowersB2Shad = [rxPowersB2Shad rxPow];
    rxPow = calcRXPower(mobilePos(index), startCenter, frequency, 'path_loss_exponent_2.9_shadowing_95_conf');
    rxPowersB1ShadWorst95 = [rxPowersB1ShadWorst95 rxPow];
    rxPow = calcRXPower(mobilePos(index), endCenter, frequency, 'path_loss_exponent_2.9_shadowing_95_conf');
    rxPowersB2ShadWorst95 = [rxPowersB2ShadWorst95 rxPow];
    
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

%% 5. Add fading effects for each point in time
% Create fading envelope
N = 65536;
fm = round(metric_speed/((3e8/frequency))); %Hz
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
time_series = 10*log10(time_series);

% Create new rx power variables
rxPowersB1Fade = [];
rxPowersB2Fade = [];
rxPowersB1ShadFade = [];
rxPowersB2ShadFade = [];
for i=1:numel(rxPowersB1)
    rxPowersB1Fade = [rxPowersB1Fade rxPowersB1(i)+...
        time_series(round((i/waveform_time)*numel(time_series)))];
end
for i=1:numel(rxPowersB2)
    rxPowersB2Fade = [rxPowersB2Fade rxPowersB2(i)+...
        time_series(round((i/waveform_time)*numel(time_series)))];
end
for i=1:numel(rxPowersB1Shad)
    rxPowersB1ShadFade = [rxPowersB1ShadFade rxPowersB1Shad(i)+...
        time_series(round((i/waveform_time)*numel(time_series)))];
end
for i=1:numel(rxPowersB2Shad)
    rxPowersB2ShadFade = [rxPowersB2ShadFade rxPowersB2Shad(i)+...
        time_series(round((i/waveform_time)*numel(time_series)))];
end



%% 6. Plot received signal power over time for all cases

figure(2);
plot(rxPowersB1, 'LineStyle', '-', 'Color', 'red');
hold on;
plot(rxPowersB2, 'LineStyle', '-', 'Color', 'blue');
title({'Case i) (No Shadowing && No Fading)', 'Received Power For Handoff'});
ylabel("Power (dBm)");
xlabel("Time (seconds)");
xlim([1, numel(rxPowersB1)]);
legend("Power from BS A", "Power from BS B");


figure(3);
plot(rxPowersB1Shad, 'LineStyle', '--', 'Color', 'red');
hold on;
plot(rxPowersB2Shad, 'LineStyle', '--', 'Color', 'blue');
title({'Case ii) (With Shadowing && No Fading)', 'Received Power For Handoff'});
ylabel("Power (dBm)");
xlabel("Time (seconds)");
xlim([1, numel(rxPowersB1)]);
legend("Power from BS A with Shadowing", "Power from BS B with Shadowing");

figure(4);
plot(rxPowersB1Fade, 'LineStyle', '-', 'Color', 'red');
hold on;
plot(rxPowersB2Fade, 'LineStyle', '-', 'Color', 'blue');
title({'Case iii) (No Shadowing && With Fading)', 'Received Power For Handoff'});
ylabel("Power (dBm)");
xlabel("Time (seconds)");
xlim([1, numel(rxPowersB1)]);
legend("Power from BS A with Fading", "Power from BS B with Fading");


figure(5);
plot(rxPowersB1ShadFade, 'LineStyle', '--', 'Color', 'red');
hold on;
plot(rxPowersB2ShadFade, 'LineStyle', '--', 'Color', 'blue');
title({'Case iv) (With Shadowing && With Fading)', 'Received Power For Handoff'});
ylabel("Power (dBm)");
xlabel("Time (seconds)");
xlim([1, numel(rxPowersB1)]);
legend("Power from BS A with Shadowing and Fading", "Power from BS B with Shadowing and Fading");


% %% 6. Determine the threshold to begin handoff
% %Find when Power from A intersects the threshold
% timeIntersection = find(rxPowersB1<lowSignalThresh, 1, 'first');
% timeIntersectionShad = find(rxPowersB1ShadWorst95<lowSignalThresh, 1, 'first');
% 
% %Subtract how long it takes for a handoff to occur
% handoffTime = 4.5;
% timeInitiateHandoff = timeIntersection - handoffTime;
% timeInitiateHandoffShad = timeIntersectionShad - handoffTime;
% 
% %Determine the threshold now
% powerAtTime = rxPowersB1(floor(timeInitiateHandoff));
% powerAtTimeShad = rxPowersB1ShadWorst95(floor(timeInitiateHandoffShad));
% 
% figure(3);
% plot(rxPowersB1, 'LineStyle', '-', 'Color', 'red');
% hold on;
% plot(rxPowersB2, 'LineStyle', '-', 'Color', 'blue');
% lowSignalThresh = ones(1, numel(rxPowersB1))*powerCallMin;
% plot(lowSignalThresh, 'LineWidth', 2.5, 'Color', 'magenta');
% handOffThresh = ones(1, numel(rxPowersB1))*powerAtTime;
% plot(handOffThresh, 'LineWidth', 2.5, 'Color', 'cyan');
% title(strcat("RXPower for Handoff \Delta_{min} = ",num2str(powerAtTime-powerCallMin), " dB"));
% ylabel("Power (dBm)");
% xlabel("Time (seconds)");
% xlim([1, numel(rxPowersB1)]);
% legend("Power from BS A", "Power from BS B",...
%     "Signal Loss Threshold", "Handoff Threshold");
% 
% figure(4);
% plot(rxPowersB1ShadWorst95, 'LineStyle', '-.', 'Color', 'red');
% hold on;
% plot(rxPowersB2ShadWorst95, 'LineStyle', '-.', 'Color', 'blue');
% lowSignalThresh = ones(1, numel(rxPowersB1))*powerCallMin;
% plot(lowSignalThresh, 'LineWidth', 2.5, 'Color', 'magenta');
% handOffThresh = ones(1, numel(rxPowersB1))*powerAtTimeShad;
% plot(handOffThresh, 'LineWidth', 2.5, 'Color', 'cyan');
% title(strcat("RXPower for Handoff : 95% Confidence Shadowing \Delta_{min} = ",num2str(powerAtTimeShad-powerCallMin), " dB"));
% ylabel("Power (dBm)");
% xlabel("Time (seconds)");
% xlim([1, numel(rxPowersB1)]);
% legend("Power from BS A Shadowing Worst Case 95% Confidence",...
%     "Power from BS A Shadowing Worst Case 95% Confidence",...
%     "Signal Loss Threshold", "Handoff Threshold");


