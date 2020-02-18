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
radius = 100; %Meters
N = iValue^2 + iValue*jValue + jValue^2;
D = sqrt(3*N)*radius;
numberOfClusters = 1;
letterCountsPerFrame = [];
sectorCountsPerLetterPerFrame = [];

baseStationLabels = {};
frequency = 2.4*10^9;
sectoringAngle = 90;
sectorZones = linspace(0, 2*pi, (360/sectoringAngle)+1);
    

%% 1.Create data structure for matlab movie
video_flag = 0; % 0 (NO) || 1 (YES)
numFrames = 65;
if video_flag
	videoName = 'many_users_sectorized';
	movieFrames = moviein(numFrames);
	frame_interf = [];
end

% Note that the coordinate system is x+j*y
%  i.e. X is real axis, and Y is imaginary axis
%  Remember to be consistent with your units

%% 2. Set location of the mobile user for each frame in the movie
%   (gives the illusion of user movement)
positionMax = 150;
mobilePos = linspace( -positionMax-positionMax*j, positionMax+10+(positionMax+10)*j, numFrames );

%For Homework 3
%Let's make a ton of mobile users!
numberOfUsers = 20;
mobilePositions = [mobilePos];
for user = 1:(numberOfUsers-1)
    ranDom = rand(1,4)*(2*positionMax)-positionMax;
    mobilePositions = [mobilePositions; linspace( ranDom(1)+ranDom(2)*j, ranDom(3)+ranDom(4)*j, numFrames )];
end


%% 3. Create the animation
for index = 1:numFrames    % Draw each frame in the movie
    clf;
    hold on;
    axis off;
    
    % Draw the serving cell and label it
%     drawCell( 0, 100, 'A_1' );
    clusterCenters = drawManyClusters(0, iValue, jValue, radius, numberOfClusters);
    % Draw the sectors! For Homework 4
    for center = clusterCenters
        for sector = sectorZones
            if abs(sin(sector)) == 1
                line( [real(center) real(center)+cos(sector)*sqrt(3)/2*radius],...
                [imag(center) imag(center)+sin(sector)*sqrt(3)/2*radius],...
                'LineStyle','--');
            else
                line( [real(center) real(center)+cos(sector)*radius],...
                [imag(center) imag(center)+sin(sector)*radius],...
                'LineStyle','--');
            end
        end
    end
    
    letterCounts = zeros(1, N);
    sectorLetterCounts = zeros(numel(sectorZones)-1, N);
    
    % Draw the mobile user(s) at the appropriate location
%     plot( mobilePos(index), 'x' );
    hold on;
    for user = 1:numberOfUsers
        %Draw all of our users
        plot( mobilePositions(user, index), 'bx' );
        [servingCellNumber, servingTierNumber, servingCellCenter] = ...
            findServingCell( mobilePositions(user, index), clusterCenters);
        line( [real(servingCellCenter) real(mobilePositions(user, index))],...
            [imag(servingCellCenter) imag(mobilePositions(user, index))],...
            'Color','blue',...
            'LineStyle','-',...
            'LineWidth', 1);
        
        %Keep track of how many users are in each base station
        letterCounts(servingTierNumber) = letterCounts(servingTierNumber)+1;
        change_imag = (imag(mobilePositions(user, index))...
            -imag(servingCellCenter));
        change_real = (real(mobilePositions(user, index))...
            -real(servingCellCenter));
        angleWithServingCell = atan(change_imag/change_real);
        if change_real <= 0
            if change_real == 0 && change_imag >= 0
                angleWithServingCell = angleWithServingCell;
            end
            angleWithServingCell = angleWithServingCell + pi;
        end
        if angleWithServingCell < 0
            angleWithServingCell = angleWithServingCell + 2*pi;
        end
        sectorNum = find(angleWithServingCell <= sectorZones, 1)-1;
        sectorLetterCounts(sectorNum, servingTierNumber) = sectorLetterCounts(sectorNum, servingTierNumber) + 1;
    end
    letterCountsPerFrame = [letterCountsPerFrame; letterCounts];
    sectorCountsPerLetterPerFrame(index ,:,:) = sectorLetterCounts;
    hold off;
    
    
    baseStationLabels(numel(baseStationLabels)+1, 1) = {convertStringsToChars(letters(servingCellNumber)+"_"+num2str(servingTierNumber))};

      
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

%% 5. Plot how many user are in each cell and sector (homework 4) over time

figure(5);
hold on;
for numCells=1:size(letterCountsPerFrame,2)
    if numCells < 4
        subplot(3, 3, numCells)
    elseif numCells == 4
        subplot(3, 3, 5)
    else
        subplot(3, 3, numCells+2)
    end
    plot(1:numFrames, letterCountsPerFrame(:,numCells)')
    title(strcat("Mobile Users in Cell ",letters(numCells)));
    ylim([0 max(max(letterCountsPerFrame))]);
    xlabel("Frame #");
    ylabel("User Count");
end
ax = gca;   % Important

figure(6);
hold on;
for numCells=1:size(letterCountsPerFrame,2)
    if numCells < 4
        subplot(3, 3, numCells)
    elseif numCells == 4
        subplot(3, 3, 5)
    else
        subplot(3, 3, numCells+2)
    end
    hold on;
    for sector=1:(numel(sectorZones)-1)
        plot(1:numFrames, sectorCountsPerLetterPerFrame(:, sector, numCells));
    end
    title(strcat("Cell ",letters(numCells), " Sector User #"));
    xlabel("Frame #");
    ylabel("User Count");
end
subplot(3,3,6)
plot(0,0,  0,0,  0,0,  0,0)
axis off
legend('Sector 1','Sector 2','Sector 3','Sector 4')