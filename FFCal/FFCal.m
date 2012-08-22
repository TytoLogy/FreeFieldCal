%------------------------------------------------------------------------
% FFCal.m
%------------------------------------------------------------------------
%
% Calibrates free-field speaker array using Bruel and Kjaer reference 
% microphone, positioned using 2-D robot controlled by RX5
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sharad.shanbhag@einstein.yu.edu
%------------------------------------------------------------------------
% Created: 2007
% (re)Created: 8 June, 2009
%
% Revisions:
%
% 	30 Oct, 07: test run: no sound, robot, etc: 433 sec. time.
% 	10 Jan, 08: Major (!) Modifications made to do tests of calibration.
%	June, 09:	Another major overhaul, now performs more tests
%	18 June, 2009 (SJS): Now makes copy of the FFCal_settings.m file in
%								the destination data directory.
%	28 October, 2009 (SJS): addeed some comments
%	18 March, 2010 (SJS):
% 		-	made changes to comply with modifications made to some of the RP
% 			toolbox functions
% 		-	updated comments/documentation
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

global DEBUG;
DEBUG = 0;

warning off MATLAB:log:logOfZero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the settings and constants for FFcal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% calibration settings
	FFCal_settings;
	% file i/o settings
	FFCal_filesettings;

	% cancel if ABORT is set
	if ABORT
		return;
	end
	
	% save a copy of the FFCal_settings.m file in the data directory
	srcfile = 'FFCal_settings.m';
	destfile = [datadir filesep srcfile];
	copyfile(srcfile, destfile);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the TDT circuits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	FFCal_tdtinit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the robot GUI
%	servo or handles.servo is the robot controller
%	interface via the TDT RX5 digital IO port
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Launch the Robot Control GUI and store the figure handle in robot_fig
	robot_fig = RoboGui;
	% retreive the handles for the robot_fig
	roboth = guihandles(robot_fig);
	handles = guidata(robot_fig);
	% setup the servo controller in handles
	handles.servo = servo;
	handles.servoinit = 1;
	% store handles
	guidata(robot_fig, handles);
	% load robot calibration data
 	load robocal
	handles = guidata(robot_fig);
	
	% Issue triggers to the servos to start and reset them
	RPtrig(servo, 1);
	RPtrig(servo, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the start and end bins for the Golay calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	start_bin = ms2bin(cal.StimDelay, indev.Fs);
	if ~start_bin
		start_bin = 1;
	end
	end_bin = ms2bin(cal.StimDelay, indev.Fs) + length(cal.A);
	zerostim = zeros(ms2bin(cal.StimDelay + cal.StimDuration, indev.Fs), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store the TDT control objects for now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	cal.indev = indev;
	cal.outdev = outdev;
	cal.servo = servo;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get a highpass filter for processing the RMS data
%  (needed to wait until Fs is known for indev)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Nyquist frequency
	cal.fnyq = indev.Fs/2;
	% filter coefficients
	[cal.coeffb, cal.coeffa] = butter(cal.forder, cal.fcutoff/cal.fnyq, 'high');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now initiate sweeps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create data directory (if it doesn't exist)
	mkdir(datadir);

	if init_spkrinfo
		% initialize controls/plot
		roboaz = 0.001 * 0.5 * interp1(azs, azdata, 0);
		set(roboth.AzimuthCtrl, 'Value', roboaz);
		set(roboth.az_value, 'String', sprintf('%.3f', roboaz));
		disp('Please press SET POSITION');
		disp('ALIGN ROBOT BASE WITH 0 AZIMUTH');
		disp('press key to continue');
		pause;
	else
		% resume at previous location
		disp('Beginning calibration where we left off.');
		disp('press key to continue');
		pause;
	end	

	% wait for user to start the acquisition
	disp('Ready to Start...');
	disp('Press key to start calibration');
	pause;

	% initialize some variables to keep track of where we are
	freq_index = 1;
	current_az = 0;
	current_el = 0;
	
	% create some data figures
	stimrespfig1 = figure;
	set(gcf, 'Name', 'Stimulus & Response 1');
	stimrespfig2 = figure;
	set(gcf, 'Name', 'Stimulus & Response 2');
	xferfig1 = figure;
	xferfig2 = figure;	
	splcurvefig = figure;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% loop through the speakers - start_speaker is defined in 
	% FFCal_settings to allow starting from a different starting speaker
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	for spkr=start_speaker:length(spkrinfo)
		speakercal.info = spkrinfo{spkr};
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% get the D/A output channel for the current speaker
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		outdev.channel = spkrinfo{spkr}.dachannel;		
		az = spkrinfo{spkr}.azimuths;
		el = spkrinfo{spkr}.elevations;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  		% Set the robot position & Set the MPX channel
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp('...moving robot...');
  		FFCal_moverobot;
		
		disp('...setting GS3 speaker channel...');
		gs3obj = gs3_open('COM6');
		status = set_speaker(gs3obj, spkrinfo{spkr}.azimuths, spkrinfo{spkr}.elevations);
		gs3_close(gs3obj);
		delete(gs3obj);
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Run the Golay (Freq/Phase) Calibration curve
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp('Now running calibration...');
		FFCal_GolayResponse;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Run the SPL curve
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		FFCal_SPLcurve;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Run the click response
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 		FFCal_ClickResponse;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Run the SPL - BW curve
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if TEST_SPL_BANDWIDTH
			% check if figures for this function exist
			% if not, create them
			if ~exist('bw_datafig', 'var')
				bw_datafig = figure;
			end
			if ~exist('bw_semilog', 'var')
				bw_semilog = figure;
			end
			if ischar(SPL_BANDWIDTH_SPEAKERS)
				if strcmp(TONE_RESPONSE_SPEAKERS, 'ALL')
					FFCal_SPLbandwidth;
				end
			elseif find( (az == SPL_BANDWIDTH_SPEAKERS(:, 1)) & (el == SPL_BANDWIDTH_SPEAKERS(:, 2)))
				FFCal_SPLbandwidth;
			end
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Run the SPL - Fctr curve
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if TEST_BW_CENTERFREQ
			% check if figures for this function exist
			% if not, create them
			if ~exist('fctrsplfig', 'var')
				fctrsplfig = figure;
			end
			if ischar(BW_CENTERFREQ_SPEAKERS)
				if strcmp(BW_CENTERFREQ_SPEAKERS, 'ALL')
					FFCal_SPLcenterfreq;
				end
			elseif find( (az == BW_CENTERFREQ_SPEAKERS(:, 1)) & (el == BW_CENTERFREQ_SPEAKERS(:, 2)))
					FFCal_SPLcenterfreq;
			end
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Run the tone response
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if TEST_TONE_RESPONSE
			% check if figures for this function exist
			% if not, create them
			if ~exist('tonerespfig', 'var')
				tonerespfig = figure;
			end
			if ischar(TONE_RESPONSE_SPEAKERS)
				if strcmp(TONE_RESPONSE_SPEAKERS, 'ALL')
					FFCal_ToneResponse;
				end
			elseif find( (az == TONE_RESPONSE_SPEAKERS(:, 1)) & (el == TONE_RESPONSE_SPEAKERS(:, 2)))
					FFCal_ToneResponse;
			end
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Run the tone-SPL response
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 		if TEST_TONE_SPL
			% check if figures for this function exist
			% if not, create them
			if ~exist('tonesplfig', 'var')
				tonesplfig = figure;
			end
			if ischar(TONE_SPL_SPEAKERS)
				if strcmp(TONE_SPL_SPEAKERS, 'ALL')
					FFCal_ToneSPLResponse;
				end
			elseif find( (az == TONE_SPL_SPEAKERS(:, 1)) & (el == TONE_SPL_SPEAKERS(:, 2)))
					FFCal_ToneSPLResponse;
			end
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Save the data for this speaker
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp('Saving data for speaker...')
		spkrinfo{spkr}.calibrated = 1;
		netcal{spkr} = caldata;
		spkrinfo{spkr}.robot_az = get(roboth.AzimuthCtrl, 'Value');
		spkrinfo{spkr}.robot_el = get(roboth.ElevationCtrl, 'Value');
		save([datadir '\' spkrinfo{spkr}.datafile], 'speakercal');
		save([datadir '\netcal-' spkrinfo{spkr}.datafile '.mat'], 'netcal');
		save spkrinfo.mat spkrinfo
		disp('...done');
	end %spkr loop end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% save the speaker calibration data in speakercal.mat
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	save speakercal.mat speakercal caldata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% need to remove the RPs to save
	cal.indev = rmfield(cal.indev, 'C');
	cal.outdev = rmfield(cal.outdev, 'C');
	cal.servo = rmfield(cal.servo, 'C');
	save([datadir '\ArrayCalData'], 'cal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exit gracefully (remove serial objects, etc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	FFCal_exit;
	warning on MATLAB:log:logOfZero

