%------------------------------------------------------------------------
% FFCal_settings.m
%------------------------------------------------------------------------
% sets up variables for FFCal
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sharad.shanbhag@einstein.yu.edu
%------------------------------------------------------------------------
% Created: 2007
% 
% Revisions:
%	10 Jan 2007 (SJS):
% 		-	moved time variables (StimInterval, StimDelay, etc) into
%			the cal.*** structure for safe keeping.
%
% 	9 June, 2009:	major revisions
%	18 March, 2010: updated comments
%--------------------------------------------------------------------------

disp('...general setup starting...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the test options settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set TEST_SPL_BANDWIDTH to 1 if you wish to run a test of 
% SPL vs. signal bandwidth.  You can also run a subset of speakers by 
% specifying the azimuth & elevation of the speakers, or for all speakers
% by setting the _SPEAKERS variable to 'ALL'
%
% for example, to test speakers at 0, 0 and 40, 0 :
%		TEST_SPL_BANDWIDTH = 1;
% 		SPL_BANDWIDTH_SPEAKERS = [ 0 0; 40 0;];
TEST_SPL_BANDWIDTH = 1;
SPL_BANDWIDTH_SPEAKERS = [ 0 0; -90 0; 0 40; 0 -40; 90 0;];

% if TEST_BW_CENTERFREQ == 1, test BW and center frequency at the 
% speakers listed in BW_CENTERFREQ_SPEAKERS
TEST_BW_CENTERFREQ = 1;
BW_CENTERFREQ_SPEAKERS = [ 0 0; -90 0; 0 40; 0 -40; 90 0;];

% if TEST_TONE_RESPONSE == 1, test tone frequency vs. amplitude at the 
% speakers listed in TONE_RESPONSE_SPEAKERS
TEST_TONE_RESPONSE = 1;
TONE_RESPONSE_SPEAKERS = 'ALL';

% if TEST_TONE_SPL == 1, test SPL vs. tone level for
% speakers listed in TONE_SPL_SPEAKERS
TEST_TONE_SPL = 1;
TONE_SPL_SPEAKERS = [ 0 0; -90 0; 0 40; 0 -40; 90 0;];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the stimulus/acquisition settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Stimulus Interval (ms)
	cal.StimInterval = 0;
	% Stimulus Duration (ms)
	cal.StimDuration = 800;
	% Duration of epoch (ms)
	cal.SweepDuration = 1000;
	% Delay of stimulus (ms)
	cal.StimDelay = 50;
	% Total time to acquire data (ms)
	cal.AcqDuration = cal.SweepDuration;
	% Total sweep time = sweep duration + inter stimulus interval (ms)
	cal.SweepPeriod = cal.SweepDuration + cal.StimInterval;
	% Stimulus ramp on/off time
	cal.StimRamp = 0;
	cal.TestRamp = 5;

	% Speaker to Calibrate (output channel on RX8(1)
	% 12-7-2008, SJS: output channel data will be determined during
	% the speaker calibration loop; this setting ignored
	cal.OutputChannel = 1;

	% Channel on RX8(2) that input is connected to
	InputChannel = 1;
	cal.InputChannel = 3;

	%Filter?
	cal.InputFilter = 1;
	%Input Filter Fc
	cal.InputFc = 100;

	%TTL pulse duration (msec)
	cal.TTLPulseDur = 1;

	% # of sweeps, usually kept at 1
	cal.NSweeps = 1;

	% # of input channels
	Nchannels = 1;
	cal.Nchannels = Nchannels;

	% Number of speakers
	Nspeakers = 144;
	cal.Nspeakers = Nspeakers;

	% lowest and highest frequency to keep in calibration data
	MINCALFREQ = 300;
	MAXCALFREQ = 15000;	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% highpass filter for processing the RMS data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% highpass cutoff frequency (Hz)
	cal.fcutoff = 200;
	% filter order
	cal.forder = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define some constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Bruel and Kjaer reference mic sensitivity (Volts/Pascal)
	BK_sense = 3.16;
	% determine the V -> Pa conversion factor
	VtoPa = 1/BK_sense;
	% precompute the volts -> RMS conversion factor for sinusoids (0.7071)
	RMSsin = 1/sqrt(2);

	% decimation factor (to shrink arrays)
	DFACT = 10;	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the calibration stimuli (Golay Code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Number of times to present Golay stim
	cal.NCalreps = 4;

	% set up the calibration iterations
	cal.Niterations = 14;
	% generate the golay pair
	[cal.A, cal.B] = golay_pair(cal.Niterations);
	% Make sure the first and last points in the output arrays are zero 
	% (D/A will remain at last voltage value in the output arrays)
	Astim = [0 cal.A 0 0];
	Bstim = [0 cal.B 0 0];
	Alen = length(cal.A);
	Blen = length(cal.B);

	% storage variables for Golay data
	speakercal.Ar = cell(cal.NCalreps, 1);
	speakercal.Br = cell(cal.NCalreps, 1);

	% RMS value for calibration test stimuli
	%	need to specify arbitrary value because the spl calibration
	%	curve data are not available yet!
	UncalibratedGolayTestRMS = 0.5;
	CalibratedGolayTestRMS = 0.5;
	Gtest_lo = 500;
	Gtest_hi = 12000;
	
	% settings for smoothing the transfer function
	SMOOVE_METHOD = 'rloess';
	SMOOVE_SPAN = 0.01;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for Vrms - dB SPL Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	cal.splminfreq = 400;
	cal.splmaxfreq = 12000;
	cal.SPLRamp = 1;

	% RMS cal reps
	cal.RMSreps = 3;

	% RMS voltage range for SPL calibration
	cal.Vrms = [(0:0.01:0.05) (0.06:0.02:0.1) (0.2:0.2:1.6)];
	cal.Nrms = length(cal.Vrms);

	% Stimulus Interval (ms)
	cal.SPLInterval = 0;
	% duration of rms SPL stimulus 
	cal.SPLDuration = 500;

	% delay for TDT RMS output
	cal.SPLDelay = 50;
	% Total time to acquire data (ms)
	cal.SPLAcquisition = cal.SPLDuration + 2*cal.SPLDelay;
	% Total sweep time = sweep duration + inter stimulus interval (ms)
	cal.SPLSweepPeriod = cal.SPLAcquisition + cal.SPLInterval;

	% storage variables for SPL-RMS data
	speakercal.SPL = cell(cal.Nrms, cal.RMSreps);

	% Desired test stimulus level (dB SPL)
	cal.SPLtestlevel_db = 40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for click testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% # of click reps to present
	cal.Nclicks = 3;

	% Stimulus Interval (ms)
	cal.ClickInterval = 0;
	% amplitude
	cal.ClickAmplitude = 3;
	% duration of click stimulus (actual click width is 1 sample)
	cal.ClickDuration = 200;

	% delay for TDT click output
	cal.ClickDelay = 50;
	% this is the delay in the click stim vector (usually 0, use TDT to set
	% actual delay)
	cal.ClickStimDelay = 0;	
	% Total time to acquire data (ms)
	cal.ClickAcquisition = cal.ClickDuration + 2*cal.ClickDelay;
	% Total sweep time = sweep duration + inter stimulus interval (ms)
	cal.ClickSweepPeriod = cal.ClickAcquisition + cal.ClickInterval;

	% storage variables for Click data
	speakercal.Click = cell(cal.Nclicks, 1);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for SPL-bandwidth testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if TEST_SPL_BANDWIDTH
		cal.bwctrfreq = 6000;
		cal.bwfreqs = [1 10 50 100 250 500 800 1000 5000 10000];
		cal.Nbwfreqs = length(cal.bwfreqs);
	
		% RMS cal reps
		cal.SPL_BWreps = 2;

		% RMS voltage range for SPL calibration
		cal.SPL_BWrms = [0.05 0.1 0.25 0.5 0.75];
		cal.Nsplbwrms = length(cal.SPL_BWrms);

		% allocate storage variables for SPL-BW-RMS data
		speakercal.SPL_BW = cell(cal.Nsplbwrms, cal.Nbwfreqs);		
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for SPL-Fctr testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if TEST_BW_CENTERFREQ || TEST_SPL_BANDWIDTH
		cal.BWbroad.fctr = 6000;
		cal.BWbroad.bw = 10000;

		cal.BWnarrow.fctr = 500:1000:10000;
		cal.BWnarrow.bw = 500;
	
		% desired test level
		cal.BWtestlevel = 40;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for tone testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if TEST_TONE_RESPONSE
		cal.ToneMinFreq = 500;
		cal.ToneMaxFreq = 9500;
		cal.ToneStep = 1000;
		cal.ToneFreqs = cal.ToneMinFreq:cal.ToneStep:cal.ToneMaxFreq;
		cal.Ntonefreqs = length(cal.ToneFreqs);

		cal.ToneLevel = 40;	% level in dB SPL

		cal.ToneDuration = 200;
		cal.ToneDelay = 25;
		cal.ToneInterval = 0;
		cal.ToneReps = 3;
		cal.ToneRamp = 5;

		% Total time to acquire data (ms)
		cal.ToneAcquisition = cal.ToneDuration + 2*cal.ToneDelay;
		% Total sweep time = sweep duration + inter stimulus interval (ms)
		cal.ToneSweepPeriod = cal.ToneAcquisition + cal.ToneInterval;

		% storage variables for Tone Data
		% array for storing magnitudes
		speakercal.ToneResp.mag = zeros(cal.Ntonefreqs, cal.ToneReps);
		% array for storing phase
		speakercal.ToneResp.phase = zeros(cal.Ntonefreqs, cal.ToneReps);
		% array for storing distortion (2nd harmonic)
		speakercal.ToneResp.distmag = zeros(cal.Ntonefreqs, cal.ToneReps);
		speakercal.ToneResp.distphase = zeros(cal.Ntonefreqs, cal.ToneReps);
		% array for storing scaling factors
		speakercal.ToneResp.scale = zeros(cal.Ntonefreqs, 1);
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for tone SPL testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if TEST_TONE_SPL
		cal.ToneSPLMinFreq = 500;
		cal.ToneSPLMaxFreq = 10000;
		cal.ToneSPLFreqStep = 500;
		cal.ToneSPLFreqs = cal.ToneSPLMinFreq:cal.ToneSPLFreqStep:cal.ToneSPLMaxFreq;
		cal.Ntonesplfreqs = length(cal.ToneSPLFreqs);

		% RMS voltage range for SPL calibration
	 	cal.ToneSPLRMS = [0.05 0.25 0.75 1.25 1.75 2.0];
		cal.Ntonesplrms = length(cal.ToneSPLRMS);

		cal.ToneSPLDuration = 200;
		cal.ToneSPLDelay = 25;
		cal.ToneSPLInterval = 0;
		cal.ToneSPLReps = 3;
		cal.ToneSPLRamp = 1;

		% Total time to acquire data (ms)
		cal.ToneSPLAcquisition = cal.ToneSPLDuration + 2*cal.ToneSPLDelay;
		% Total sweep time = sweep duration + inter stimulus interval (ms)
		cal.ToneSPLSweepPeriod = cal.ToneSPLAcquisition + cal.ToneSPLInterval;

		% storage variables for Tone SPL Data
		% array for storing magnitudes
		speakercal.ToneSPLData = cell(cal.Ntonesplrms, cal.Ntonesplfreqs, cal.ToneSPLReps);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDT Circuit Filenames and paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Tytology root directory
	tytologyroot = 'H:\Code\TytoLogy\';
	
	% output device configuration
	outdev.Circuit_Path = [tytologyroot 'toolbox\TDT\Circuits\RX8\50KHz\'];
% 	outdev.Circuit_Name = 'RX8_1_SingleChannelOutput';
	outdev.Circuit_Name = 'RX8_Array_SingleChannelOutput';
	% Dnum = device number - this is for RX8 (1)
	outdev.Dnum=1;
	outdev.C = [];
	outdev.status = 0;

	% input device configuration
	if cal.InputFilter
		indev.Circuit_Path = [tytologyroot 'Toolbox\TDT\Circuits\RX8_2\50KHz\'];
		indev.Circuit_Name = 'RX8_2_SingleChannelInput_HPFilt';	
	else
		indev.Circuit_Path = [tytologyroot 'Toolbox\TDT\Circuits\RX8_2\50KHz\'];
		indev.Circuit_Name = 'RX8_2_SingleChannelInput';
	end
	% Dnum = device number - this is for RX8 (2)
	indev.Dnum=2;
	indev.C = [];
	indev.status = 0;
	
	% servo (robot) controller device configuration
	servo.Circuit_Path = [tytologyroot 'Toolbox\TDT\Circuits\RX5\'];
	servo.Circuit_Name = 'ServoPulse';
	% Dnum = device number - this is for RX5 (1)
	indev.Dnum=2;
	servo.C = [];
	servo.status = 0;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check to make sure the circuits exist!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tmpfile = [fullfile(outdev.Circuit_Path, outdev.Circuit_Name) '.rco'];
	if ~exist(tmpfile, 'file')
		error([mfilename ': output circuit ' tmpfile ' not found']);
	end
	tmpfile = [fullfile(indev.Circuit_Path, indev.Circuit_Name) '.rco'];
	if ~exist(tmpfile, 'file')
		error([mfilename ': input circuit ' tmpfile ' not found']);
	end
	tmpfile = [fullfile(servo.Circuit_Path, servo.Circuit_Name) '.rco'];
	if ~exist(tmpfile, 'file')
		error([mfilename ': servo circuit ' tmpfile ' not found']);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify i/o channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	outdev.channel = cal.OutputChannel;
	indev.channel = cal.InputChannel;	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plot labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

