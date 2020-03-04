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
radius = (2*(1000))/sqrt(3); %Meters
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
    

%% 1.Create data structure for matlab movie
video_flag = 1; % 0 (NO) || 1 (YES)
distance = 2000; % Meters
speed = 50; % Miles per hour
metric_speed = speed / 2.237; % Meters per second
% If each frame represents a second
numFrames = 20;
if video_flag
	videoName = 'steering_video';
	movieFrames = moviein(numFrames);
	frame_interf = [];
end
if video_flag
	videoName2 = 'steering_DOA_video';
	movieFrames2 = moviein(numFrames);
	frame_interf2 = [];
end

%% 2. Set location of the mobile user for each frame in the movie
%   (gives the illusion of user movement)
positionMax = 150;
startCenter = (1.0e+03*(-0.5 + 1.7*j))/2;
endCenter = (1.0e+03*(.834 - .2*i))/2;
mobilePos = linspace(endCenter,0-startCenter, numFrames );
% mobilePos = ones(1, numFrames)*(1e2)*(3.7834 - 2*i);

%For Homework 3
%Let's simulate handoff!


%% 3. Create the animation
for index = 1:numFrames    % Draw each frame in the movie
    clf;
    hold on;
    axis off;
    
    % Draw the serving cell and label it
    drawCell( 0, radius, 'A' );
    clusterCenters = 0;
    axis equal;
%     clusterCenters = drawManyClusters(0, iValue, jValue, radius, numberOfClusters);
    letterCounts = zeros(1, N);
    
    % Draw our antenna locations
%     for i=1:num_antennas
%        plot(antenna_locations(i), 'g^'); 
%     end
    
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
    rxPow = calcRXPower(mobilePos(index), startCenter, frequency, 'path_loss_exponent_2.9');
    rxPowersB1 = [rxPowersB1 rxPow];
    rxPow = calcRXPower(mobilePos(index), startCenter, frequency, 'path_loss_exponent_2.9_shadowing');
    rxPowersB1Shad = [rxPowersB1Shad rxPow];
    rxPow = calcRXPower(mobilePos(index), startCenter, frequency, 'path_loss_exponent_2.9_shadowing_95_conf');
    rxPowersB1ShadWorst95 = [rxPowersB1ShadWorst95 rxPow];    
    
    % Each row is a frame, each column is the corresponding antenna
%     ant_strengths = [ant_strengths; sig_x];
    
    % Capture the frame for the movie
	movieFrames(:, index) = getframe(gcf);	
    
    	
end
clf;
hold on;
axis off;

% Draw the serving cell and label it
drawCell( 0, radius, 'A' );
clusterCenters = 0;
axis equal;
%     clusterCenters = drawManyClusters(0, iValue, jValue, radius, numberOfClusters);
letterCounts = zeros(1, N);

% Draw our antenna locations
%     for i=1:num_antennas
%        plot(antenna_locations(i), 'g^'); 
%     end

% Draw the mobile user(s) at the appropriate location
plot( mobilePos(1), 'rx' );

 % Draw a line connecting the center (basestation) of the serving cell 
%    and the mobile user
line( [real(mobilePos(1)) real(servingCellCenter)],...
    [imag(mobilePos(1)) imag(servingCellCenter)],...
        'Color','blue',...
        'LineStyle','-',...
        'LineWidth', 1);
hold off;
%% 4. DOA
figure(2)
for index = 1:numFrames    % Draw each frame in the movie
    [steering, spectrum] = ulaSteering(mobilePos(index), clusterCenters, frequency, num_antennas, ant_space);
    
    fplot(spectrum/5, [0 pi]);
    title(strcat("Spatial Spectrum Frame ", num2str(index)));
    xlabel("DOA (radians)");
    ylabel("Spectrum Value");
    
    actual_DOA = atan((real(mobilePos(index))...
            -real(0))/(imag(0)...
            -imag(mobilePos(index))));
    xline(actual_DOA, 'g--');
    movieFrames2(:, index) = getframe(gcf);
end
[steering, spectrum] = ulaSteering(mobilePos(1), clusterCenters, frequency, num_antennas, ant_space);

fplot(spectrum/5, [0 pi]);
title(strcat("Spatial Spectrum Frame ", num2str(1)));
xlabel("DOA (radians)");
ylabel("Spectrum Value");

actual_DOA = atan((real(mobilePos(1))...
        -real(0))/(imag(0)...
        -imag(mobilePos(1))));
xline(actual_DOA, 'g--');

%% 5. DOA for multiple spacings

spectrums = [];
count = 0;
figure(3);
for d=linspace(0.1, 1, 10)
    count = count + 1;
    [next_steer, next_spectrum] = ulaSteering(mobilePos(1), clusterCenters, frequency, num_antennas, d);
    spectrums = [spectrums next_spectrum];
    subplot(2, 5, count)
    fplot(next_spectrum/5, [0 pi]);
    actual_DOA = atan((real(mobilePos(1))...
            -real(0))/(imag(0)...
            -imag(mobilePos(1))));
    xline(actual_DOA, 'g--');
    title(strcat("Spatial Spectrum ", num2str(d), "{\lambda}"));
    xlabel("DOA (radians)");
    ylabel("Spectrum Value");
end



%% 6. Write the movie to a file (videoName.avi)
if video_flag
	v = VideoWriter(videoName);
	% Modify as needed to make your animation faster/slower
	v.FrameRate = 10;
	open(v);
	writeVideo(v,movieFrames);
	close(v);
end

if video_flag
	v2 = VideoWriter(videoName2);
	% Modify as needed to make your animation faster/slower
	v2.FrameRate = 10;
	open(v2);
	writeVideo(v2,movieFrames2);
	close(v2);
end


