% main.m
%     Cameron J. Calv
%     Add clustering and method of finding serving cell for mobile user 
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
radius = 100; %Meters
N = iValue^2 + iValue*jValue + jValue^2;
D = sqrt(3*N)*radius;
numberOfClusters = 1;

%For tracking carrier strength
rxPowersN3 = [];
rxPowersN4 = [];
rxPowersN3Shadowing = [];
rxPowersN4Shadowing = [];

%For tracking interference strength
interfereColors = {'green', 'magenta', 'red'};
interfereLineTypes = {'-', '--'};
interferePowersN3 = [];
interferePowersN4 = [];
interferePowersN3Shadowing = [];
interferePowersN4Shadowing = [];

%Assumptions for SIR
assumpOneSIRN3 = [];
assumpOneSIRN4 = [];
assumpTwoSIRN3 = [];
assumpTwoSIRN4 = [];

baseStationLabels = {};
frequency = 2.4*10^9;
    

%% 1.Create data structure for matlab movie
video_flag = 0; % 0 (NO) || 1 (YES)
numFrames = 50;
if video_flag
	videoName = '1st_tier_interferers';
	movieFrames = moviein(numFrames);
	frame_interf = [];
end

% Note that the coordinate system is x+j*y
%  i.e. X is real axis, and Y is imaginary axis
%  Remember to be consistent with your units

%% 2. Set location of the mobile user for each frame in the movie
%   (gives the illusion of user movement)
mobilePos = linspace( -150-150*j, 150+160*j, numFrames );

%% 3. Create the animation
for index = 1:numFrames    % Draw each frame in the movie
    clf;
    hold on;
    axis off;
    
    % Draw the serving cell and label it
%     drawCell( 0, 100, 'A_1' );
    clusterCenters = drawManyClusters(0, iValue, jValue, radius, N);
    
    % Draw the mobile user at the appropriate location
    plot( mobilePos(index), 'x' );
    
    % Draw a line connecting the center (basestation) of the serving cell 
    %    and the mobile user
    [servingCellNumber, servingTierNumber, servingCellCenter] = findServingCell( mobilePos(index), clusterCenters);
    line( [real(servingCellCenter) real(mobilePos(index))],...
        [imag(servingCellCenter) imag(mobilePos(index))],...
            'Color','blue',...
            'LineStyle','-',...
            'LineWidth', 1);
    hold off;
    
    baseStationLabels(numel(baseStationLabels)+1, 1) = {convertStringsToChars(letters(servingCellNumber)+"_"+num2str(servingTierNumber))};

    
    % Determine receive power
%     rxPow = calcRXPower(mobilePos(index), servingCellCenter, frequency, 'friis';
    rxPow = calcRXPower(mobilePos(index), servingCellCenter, frequency, 'path_loss_exponent_3');
    rxPowersN3 = [rxPowersN3 rxPow];
    rxPow = calcRXPower(mobilePos(index), servingCellCenter, frequency, 'path_loss_exponent_4');
    rxPowersN4 = [rxPowersN4 rxPow];
    rxPow = calcRXPower(mobilePos(index), servingCellCenter, frequency, 'path_loss_exponent_3_shadowing');
    rxPowersN3Shadowing = [rxPowersN3Shadowing rxPow];
    rxPow = calcRXPower(mobilePos(index), servingCellCenter, frequency, 'path_loss_exponent_4_shadowing');
    rxPowersN4Shadowing = [rxPowersN4Shadowing rxPow];
    
    %Determine interference power from 1st tier co-channel cells
    interferingCellCenters = [];
    interferePowN3 = 0;
    interferePowN4 = 0;
    interferePowN3Shadowing = 0;
    interferePowN4Shadowing = 0;
    for c=1:6 %Every direction!
        %Find the center of the next interfering cell
        direction = pi/6 + (c-1)*(pi/3);
        nextPosition = servingCellCenter + 2*(sqrt(3)/2)*radius*iValue*(cos(direction)+j*sin(direction));
        direction = direction + pi/3;
        interferingCellCenter = nextPosition + 2*(sqrt(3)/2)*radius*jValue*(cos(direction)+j*sin(direction));
        
        %Calculate the path loss for each scenario
        interferePow = calcRXPower(mobilePos(index), interferingCellCenter, frequency, 'path_loss_exponent_3');
        interferePowN3 = interferePowN3 + dbm2mw(interferePow);
        interferePow = calcRXPower(mobilePos(index), interferingCellCenter, frequency, 'path_loss_exponent_4');
        interferePowN4 = interferePowN4 + dbm2mw(interferePow);
        interferePow = calcRXPower(mobilePos(index), interferingCellCenter, frequency, 'path_loss_exponent_3_shadowing');
        interferePowN3Shadowing = interferePowN3Shadowing + dbm2mw(interferePow);
        interferePow = calcRXPower(mobilePos(index), interferingCellCenter, frequency, 'path_loss_exponent_4_shadowing');
        interferePowN4Shadowing = interferePowN4Shadowing + dbm2mw(interferePow);
        
        %Draw another line from the user to the interfering cell
        line( [real(interferingCellCenter) real(mobilePos(index))],...
            [imag(interferingCellCenter) imag(mobilePos(index))],...
            'Color',interfereColors{mod(c-1,3)+1},...
            'LineStyle',interfereLineTypes{mod(c-1, 2)+1},...
            'LineWidth', 0.5);
        
    end
    interferePowersN3 = [interferePowersN3 mw2dbm(interferePowN3)];
    interferePowersN4 = [interferePowersN4 mw2dbm(interferePowN4)];
    interferePowersN3Shadowing = [interferePowersN3Shadowing mw2dbm(interferePowN3Shadowing)];
    interferePowersN4Shadowing = [interferePowersN4Shadowing mw2dbm(interferePowN4Shadowing)];
    
    % Under a few assumptions, calculate the SIR
    %%Assumption 1: Distance between all interferers is about the same
    beta = 3;
    assumpOneSIRN3 = [assumpOneSIRN3 radius^(-beta)/(6*(D^(-beta)))];
    beta = 4;
    assumpOneSIRN4 = [assumpOneSIRN4 radius^(-beta)/(6*(D^(-beta)))];

    
    %%Assumption 2: Farthest interferes (D+R distance), closest (D-R
    %%distance), intermediate (D distance)
    beta = 3;
    assumpTwoSIRN3 = [assumpTwoSIRN3 radius^(-beta)/...
        (2*((D-radius)^(-beta))+2*((D)^(-beta))+2*((D+radius)^(-beta)))];
    beta = 4;
    assumpTwoSIRN4 = [assumpTwoSIRN4 radius^(-beta)/...
        (2*((D-radius)^(-beta))+2*((D)^(-beta))+2*((D+radius)^(-beta)))];
    
    % Capture the frame for the movie
	movieFrames(:, index) = getframe(gcf);	
end

%% 4. Write the movie to a file (videoName.avi)
if video_flag
	v = VideoWriter(videoName);
	% Modify as needed to make your animation faster/slower
	v.FrameRate = 5;
	open(v);
	writeVideo(v,movieFrames);
	close(v);
end

%% 5. Plot received signal power over time 

figure(2);
plot(rxPowersN3);
hold on;
plot(rxPowersN4);
plot(rxPowersN3Shadowing);
plot(rxPowersN4Shadowing);
title("Received Signal Power");
ylabel("Power (dBm)");
xlabel("Serving Base Station");
legend("Exponent = 3 no Shadowing", "Exponent = 4 no Shadowing", "Exponent = 3 With Shadowing", "Exponent = 4 With Shadowing");
set(gca,'xtick',[1:numel(baseStationLabels)],'xticklabel',baseStationLabels)

%% 6. Plot 1st Tier Interference power over time

figure(3);
plot(interferePowersN3);
hold on;
plot(interferePowersN4);
plot(interferePowersN3Shadowing);
plot(interferePowersN4Shadowing);
title("1st-Tier Interference Signal Power");
ylabel("Power (dBm)");
xlabel("Serving Base Station");
legend("Exponent = 3 no Shadowing", "Exponent = 4 no Shadowing", "Exponent = 3 With Shadowing", "Exponent = 4 With Shadowing");
set(gca,'xtick',[1:numel(baseStationLabels)],'xticklabel',baseStationLabels)


%% 7. Plot Signal to Interference Ratio over time.

figure(4);
plot(rxPowersN3-interferePowersN3);
hold on;
plot(rxPowersN4-interferePowersN4);
plot(rxPowersN3Shadowing-interferePowersN3Shadowing);
plot(rxPowersN4Shadowing-interferePowersN4Shadowing);
plot(assumpOneSIRN3, '--')
plot(assumpOneSIRN4, '--')
plot(assumpTwoSIRN3, '--')
plot(assumpTwoSIRN4, '--')
title("Signal to Interference Ratio");
ylabel("SIR");
xlabel("Serving Base Station");
legend("Exponent = 3 no Shadowing", "Exponent = 4 no Shadowing",...
    "Exponent = 3 With Shadowing", "Exponent = 4 With Shadowing",...
    "Exponent = 3 1st Assumption", "Exponent = 4 1st Assumption",...
    "Exponent = 3 2nd Assumption", "Exponent = 4 2nd Assumption");
set(gca,'xtick',[1:numel(baseStationLabels)],'xticklabel',baseStationLabels)