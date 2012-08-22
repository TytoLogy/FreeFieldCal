% SpeakerMap
% 
% script to display information about speaker array
%
% Important:
% 	Make sure xlsmap.mat is in the path or in the 
% 	same directory as SpeakerMap.m
%
%--------------------------------------------------------------------------
% organization of xlsmap is:
% [144 X 7]
% 
% xlsmap(:, 1)	speaker id # (1-144)
% xlsmap(:, 2)	speaker azimuth (degrees, polar coordinates)
% xlsmap(:, 3)	speaker elevation (degrees, polar coordinates)
% xlsmap(:, 4)	GS3 (multiplexer) block #
% xlsmap(:, 5)	wiring terminal block channel (1,3,5,...,49), 6 blocks total
% xlsmap(:, 6)	RX8(2) D/A channel/GS3 input channel (1-24)
% xlsmap(:, 7)	GS3 output bank channel # (1-6, one for each input channel)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbha@aecom.yu.edu
%--------------------------------------------------------------------------
% Revision History
%	4 Mar 2009 (SJS): Created
%--------------------------------------------------------------------------
% To Do:
%--------------------------------------------------------------------------

clear all
close all

% load xlsmap
load xlsmap

% pull off the columns into individual vectors
% (more for readibility than for function)
SpeakerID = xlsmap(:, 1);
Azimuths = xlsmap(:, 2);
Elevations = xlsmap(:, 3);
GS3Block = xlsmap(:, 4);
Terminal = xlsmap(:, 5);
DAChannel = xlsmap(:, 6);
GS3Channel = xlsmap(:, 7);
Nspeakers = length(xlsmap);

%--------------------------------------------------------------------------
% first, make a plot of the speaker numbers as a function of location
%--------------------------------------------------------------------------

% Create figure 1 and axes
F1 = figure(1);

% get the current axis handle
A1 = gca;

% Title & Labels
title({'Free-Field Array';'Speaker Numbers and Locations';'   '})
xlabel('Azimuth (degrees)')
ylabel('Elevation (degrees)')

% set x and y axis limits
xlim([min(Azimuths)-20 max(Azimuths)]+10);
ylim([min(Elevations)-20 max(Elevations)]+10);

% set up the ticks to have increments of 10 degrees
xticks = min(Azimuths):10:max(Azimuths);
yticks = min(Elevations):10:max(Elevations);
set(A1, 'XTick', xticks);
set(A1, 'YTick', yticks);

% turn on grid
grid on;
set(A1, 'GridLineStyle', '-')
set(A1, 'LineWidth', 0.1)

% Plot the speaker numbers on the axes

% create an array to hold the text handles (will need handles to set
% properties later)
textHandles1= zeros(Nspeakers, 1);

% now loop through speakers
for n = 1:Nspeakers
	% draw the text on the plot and store handle
	textHandles1(n) = text(Azimuths(n), Elevations(n), num2str(SpeakerID(n)));

	% set some properties for better appearance
	set(textHandles1(n), 'Color', 'blue');
	set(textHandles1(n), 'FontName', 'Verdana');
	set(textHandles1(n), 'FontWeight', 'Bold');
	set(textHandles1(n), 'FontSize', 8)
	set(textHandles1(n), 'HorizontalAlignment', 'center')
end	



%--------------------------------------------------------------------------
% now, plot the speaker D/A channel numbers as a function of location
%--------------------------------------------------------------------------

% Create figure 1 and axes
F2 = figure(2);

% get the current axis handle
A2 = gca;

% Title & Labels
title({'Free-Field Array';'D-A Channels and Locations';'   '})
xlabel('Azimuth (degrees)')
ylabel('Elevation (degrees)')

% set x and y axis limits
xlim([min(Azimuths)-20 max(Azimuths)]+10);
ylim([min(Elevations)-20 max(Elevations)]+10);

% set up the ticks to have increments of 10 degrees
xticks = min(Azimuths):10:max(Azimuths);
yticks = min(Elevations):10:max(Elevations);
set(A2, 'XTick', xticks);
set(A2, 'YTick', yticks);

% turn on grid
grid on;
set(A2, 'GridLineStyle', '-')
set(A2, 'LineWidth', 0.1)

% Plot the speaker numbers on the axes

% create an array to hold the text handles (will need handles to set
% properties later)
textHandles2= zeros(Nspeakers, 1);

% now loop through speakers
for n = 1:Nspeakers
	% draw the text on the plot and store handle
	textHandles2(n) = text(Azimuths(n), Elevations(n), num2str(DAChannel(n)));

	% set some properties for better appearance
	set(textHandles2(n), 'Color', [0.75 0 0]);
	set(textHandles2(n), 'FontName', 'Verdana');
	set(textHandles2(n), 'FontWeight', 'Bold');
	set(textHandles2(n), 'FontSize', 10)
	set(textHandles2(n), 'HorizontalAlignment', 'center')
end	
