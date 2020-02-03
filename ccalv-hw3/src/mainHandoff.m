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

%% 5. Plot received signal power over time 

figure(2);
plot(rxPowersB1, 'LineStyle', '-', 'Color', 'red');
hold on;
plot(rxPowersB2, 'LineStyle', '-', 'Color', 'blue');
plot(rxPowersB1Shad, 'LineStyle', '--', 'Color', 'red');
plot(rxPowersB2Shad, 'LineStyle', '--', 'Color', 'blue');
plot(rxPowersB1ShadWorst95, 'LineStyle', '-.', 'Color', 'red');
plot(rxPowersB2ShadWorst95, 'LineStyle', '-.', 'Color', 'blue');
lowSignalThresh = ones(1, numel(rxPowersB1))*powerCallMin;
plot(lowSignalThresh, 'LineWidth', 2.5, 'Color', 'magenta');
title("Received Signal Power For Handoff");
ylabel("Power (dBm)");
xlabel("Time (seconds)");
xlim([1, numel(rxPowersB1)]);
legend("Power from BS A", "Power from BS B",...
    "Power from BS A with Shadowing", "Power from BS B with Shadowing",...
    "Power from BS A Shadowing Worst Case 95% Confidence",...
    "Power from BS A Shadowing Worst Case 95% Confidence",...
    "Signal Loss Threshold");

%% 6. Determine the threshold to begin handoff
%Find when Power from A intersects the threshold
timeIntersection = find(rxPowersB1<lowSignalThresh, 1, 'first');
timeIntersectionShad = find(rxPowersB1ShadWorst95<lowSignalThresh, 1, 'first');

%Subtract how long it takes for a handoff to occur
handoffTime = 4.5;
timeInitiateHandoff = timeIntersection - handoffTime;
timeInitiateHandoffShad = timeIntersectionShad - handoffTime;

%Determine the threshold now
powerAtTime = rxPowersB1(floor(timeInitiateHandoff));
powerAtTimeShad = rxPowersB1ShadWorst95(floor(timeInitiateHandoffShad));

figure(3);
plot(rxPowersB1, 'LineStyle', '-', 'Color', 'red');
hold on;
plot(rxPowersB2, 'LineStyle', '-', 'Color', 'blue');
lowSignalThresh = ones(1, numel(rxPowersB1))*powerCallMin;
plot(lowSignalThresh, 'LineWidth', 2.5, 'Color', 'magenta');
handOffThresh = ones(1, numel(rxPowersB1))*powerAtTime;
plot(handOffThresh, 'LineWidth', 2.5, 'Color', 'cyan');
title(strcat("RXPower for Handoff \Delta_{min} = ",num2str(powerAtTime-powerCallMin), " dB"));
ylabel("Power (dBm)");
xlabel("Time (seconds)");
xlim([1, numel(rxPowersB1)]);
legend("Power from BS A", "Power from BS B",...
    "Signal Loss Threshold", "Handoff Threshold");

figure(4);
plot(rxPowersB1ShadWorst95, 'LineStyle', '-.', 'Color', 'red');
hold on;
plot(rxPowersB2ShadWorst95, 'LineStyle', '-.', 'Color', 'blue');
lowSignalThresh = ones(1, numel(rxPowersB1))*powerCallMin;
plot(lowSignalThresh, 'LineWidth', 2.5, 'Color', 'magenta');
handOffThresh = ones(1, numel(rxPowersB1))*powerAtTimeShad;
plot(handOffThresh, 'LineWidth', 2.5, 'Color', 'cyan');
title(strcat("RXPower for Handoff : 95% Confidence Shadowing \Delta_{min} = ",num2str(powerAtTimeShad-powerCallMin), " dB"));
ylabel("Power (dBm)");
xlabel("Time (seconds)");
xlim([1, numel(rxPowersB1)]);
legend("Power from BS A Shadowing Worst Case 95% Confidence",...
    "Power from BS A Shadowing Worst Case 95% Confidence",...
    "Signal Loss Threshold", "Handoff Threshold");


