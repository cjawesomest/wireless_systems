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
N = iValue^2 + iValue*jValue + jValue^2;
numberOfClusters = 1;

rxPowers = [];
baseStationLabels = {};
frequency = 2.4*10^9;
    

%% 1.Create data structure for matlab movie
video_flag = 0; % 0 (NO) || 1 (YES)
numFrames = 20;
if video_flag
	videoName = 'homework_1_movie';
	movieFrames = moviein(numFrames);
	frame_interf = [];
end

% Note that the coordinate system is x+j*y
%  i.e. X is real axis, and Y is imaginary axis
%  Remember to be consistent with your units

%% 2. Set location of the mobile user for each frame in the movie
%   (gives the illusion of user movement)
mobilePos = linspace( -100-100*j, 100+150*j, numFrames );

%% 3. Create the animation
for index = 1:numFrames    % Draw each frame in the movie
    clf;
    hold on;
    axis off;
    
    % Draw the serving cell and label it
%     drawCell( 0, 100, 'A_1' );
    clusterCenters = drawManyClusters(0, iValue, jValue, 100, 7);
    
    % Draw the mobile user at the appropriate location
    plot( mobilePos(index), 'x' );
    
    % Draw a line connecting the center (basestation) of the serving cell 
    %    and the mobile user
    [servingCellNumber, servingTierNumber, servingCellCenter] = findServingCell( mobilePos(index), clusterCenters);
    line( [real(servingCellCenter) real(mobilePos(index))], [imag(servingCellCenter) imag(mobilePos(index))] );
    hold off;
    
    baseStationLabels(numel(baseStationLabels)+1, 1) = {convertStringsToChars(letters(servingCellNumber)+"_"+num2str(servingTierNumber))};

    
    % Determine receive power
    rxPow = calcRXPower(mobilePos(index), servingCellCenter, frequency, @friisFreeSpace);
    rxPowers = [rxPowers rxPow];
    
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

%% 5. Plot received signal power over time (Friis Free Space)

figure(2);
plot(rxPowers);
title("Received Signal Power");
ylabel("Power (watts)");
xlabel("Serving Base Station");
set(gca,'xtick',[1:numel(baseStationLabels)],'xticklabel',baseStationLabels)
