%------------------------------------------------------------------------
% BuildNetcal.m
%	Script to process free-field calibration data from FFCal program
%
%	Usage:
% 
% 		1) Collect data using FFCal
% 		2)	Set paths to data 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created:	?
% 
% Revisions:
%	10 June, 2009 (SJS): added code to compute dbspl data (regression)
%						 migrating this to a function, processFFCaldata.m
%	18 June, 2009 (SJS): renamed to BuildNetcal.m
%	27 October, 2009 (SJS): adding code for tone calibration processing
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Some setup
%------------------------------------------------------------------------
	clear all
	close all

	% turn of log of Zero warning message - need to do this to avoid being
	% inundated with warning messages
	warning off MATLAB:log:logOfZero

%------------------------------------------------------------------------
% some command "settings"
%------------------------------------------------------------------------
	% SAVE_SPEAKERDATA used to tell program to save the individual
	% speaker calibration data in an individual file
	SAVE_SPEAKERDATA = 0;
	
	% SAVE_NETDATA tells program to save the net calibration file (usually
	% used, if you want BuildNetcal to do anything useful...)
	SAVE_NETDATA = 1;

	% option and settings to smooth calibration data files
	SMOOVE = 1;
	SMOOVE_METHOD = 'rloess';
	SMOOVE_SPAN = 0.1;

%------------------------------------------------------------------------
% paths and information
%------------------------------------------------------------------------
	% this reads in the settings from the data acquisition - set it to point
	% to the appropriate file in the calibration data directory, e.g.:
	datapath = uigetdir(pwd, 'Select raw data location (FFCal output data)');
	if isequal(datapath, 0)
		disp('cancelled by user...')
		return
	end
	% build ffcalsettings file name (FFCal_settings.m)
	ffcalsettingsfile = fullfile(datapath, 'FFCal_settings.m');
	
	% look for the FFCal_settings.m file
	if exist(ffcalsettingsfile, 'file')
		% if it's there, run the script
		run(ffcalsettingsfile);
	else
		warning('%s: FFCal_settings.m file not found in directory %s', ...
						mfilename, datapath);
		[tmpfile, tmppath] = uigetfile('FFCal_settings.m', 'Open FFCal_settings.m file');
		if isequal(tmpfile, 0)
			disp('cancelled')
			return
		else
			ffcalsettingsfile = fullfile(tmppath, tmpfile);
			run(ffcalsettingsfile);
		end
	end
	
	
	% specify the input directory and input data name
	dataname = query_userstring('Enter data name (e.g., 18-Mar-2020)')
		
	% specify the output directory
	outdatapath = datapath;
	outdatapath = uigetdir(outdatapath, 'Select destination directory for processed data')
	if ~exist(outdatapath, 'dir')
		mkdir(outdatapath);
	end

	% specify output (destination) net calibration file
	netcalfile = [outdatapath filesep 'netcal-' dataname '.mat'];

	% need to load ArrayCalData.mat file generated by the 
	arraycalfile = [datapath filesep 'ArrayCalData.mat'];
	if ~exist(arraycalfile, 'file')
		warning('%s: ArrayCalData.mat file not found in directory %s', ...
						mfilename, datapath);
		[tmpfile, tmppath] = uigetfile('ArrayCalData.mat', 'Open ArrayCalData.mat file');
		if isequal(tmpfile, 0)
			disp('cancelled')
			return
		else
			arraycalfile = fullfile(tmppath, tmpfile);
		end
	end	
	load(arraycalfile);

	indev = cal.indev;
	outdev = cal.outdev;

%------------------------------------------------------------------------
% allocate variables
%------------------------------------------------------------------------
netcal = cell(Nspeakers, 1);

%------------------------------------------------------------------------
% PROCESS THE DATA!!!!
%------------------------------------------------------------------------

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% loop throught the speakers (usually 144)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for spkrIndex = 1:cal.Nspeakers
% 	for spkrIndex = 73:73
		% get the information for the current speaker
		currentspkr = speaker_data('speakernum', spkrIndex)
		az = currentspkr.azimuths
		el = currentspkr.elevations
		
		% build the speaker data file names
		fname = sprintf('%d_%d', az, el);
		% need to replace '-' signs with 'n'
		fname(find(fname=='-')) = 'n';		
		speakercalfile = fullfile(datapath, [fname '.mat'])
		speakernetfile = fullfile(outdatapath, ['netcal-' fname '.mat']);

		% load the data for this speaker
		load(speakercalfile);

		% set the caldata struct's speaker information
		caldata.speakerinfo = currentspkr;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Now, process the Golay Calibration Data
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp(sprintf('Computing Transfer function, speaker az %d el %d ...', az, el));

		for n=1:cal.NCalreps
			Ar=speakercal.Ar{n};
			Br=speakercal.Br{n};

			Arlen = length(Ar);	Brlen = length(Br);

			% Run the golay pair and the response through golay_impulse to get
			% the transfer function
			[H{n}, I{n}, Freq] = golay_impulse(cal.A, cal.B, Ar, Br, indev.Fs);

			% keep the first M points in the fft arrays (the rest are
			% duplicates)
			M = length(Freq);
			Tg{n} = normalize(abs(H{n}(1:M)));
			Pg{n} = unwrap(angle(H{n}(1:M)));

			Tgm(n, :) = Tg{n};
			Pgm(n, :) = Pg{n};
			% catch DC == 0
			if Tgm(n, 1) == 0
				Tgm(n, 1) = min(Tg{n});
			end
		end

		% take means of variables
		Tgm=mean(Tgm);
		Pgm=mean(Pgm);

		figure(1)
		% plot mean response and trials
		c=['r' 'm' 'g' 'k' 'y' 'b'];
		subplot(4, 1, 1),	plot(Freq, db(Tg{n}), c(n));
		subplot(4, 1, 2),	plot(Freq, Pg{n}, c(n));
		for n=1:cal.NCalreps
			subplot(4, 1, 1)
			hold on 
			plot(Freq, db(Tg{n}), c(n)) 
			hold off
			subplot(4, 1, 2)
			hold on
			plot(Freq, Pg{n}, c(n));
			hold off
		end
		subplot(4, 1, 1)
		hold on
		plot(Freq, db(Tgm), 'k') 
		hold off
		title(sprintf('Transfer Function Az: %d El %d', az, el));
		ylabel('Gain (dB)');
		subplot(4, 1, 2)
		hold on
		plot(Freq, Pgm, 'k');
		hold off
		xlabel('Frequency (Hz)');
		ylabel('Phase');

		% Process the transfer function data for calibration duties
		disp('...Processing net calibration data...')
		minfreq = max(find(Freq <= MINCALFREQ));
		maxfreq = min(find(Freq >= MAXCALFREQ));
		caldata.minfreq = minfreq;
		caldata.maxfreq = maxfreq;
		caldata.freq = downsample(Freq(minfreq:maxfreq), DFACT);
		caldata.mag = decimate(Tgm(minfreq:maxfreq), DFACT);
		caldata.phase = decimate(Pgm(minfreq:maxfreq), DFACT);
		magdb = db(caldata.mag);
		mindb = min(magdb);
		caldata.maginv = power(10, (mindb - magdb) ./ 20);
		if SMOOVE
			caldata.maginv = smooth(caldata.maginv, SMOOVE_SPAN, SMOOVE_METHOD)';
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Analyze the RMS-SPL  Data  for Noise
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% loop throught Nrms and RMSrepetitions
		for n = 1:cal.Nrms
			for r = 1:cal.RMSreps
				% filter the raw data, convert to Pascal, RMS and dB SPL
				splfilt = filter(cal.coeffb, cal.coeffa, speakercal.SPL{n, r});
				unfilt_rms(r) = rms(VtoPa * speakercal.SPL{n, r});
				unfilt_dbspl(r) = dbspl(unfilt_rms(r));
				filt_rms(r) = rms(VtoPa * splfilt);
				filt_dbspl(r) = dbspl(filt_rms(r));
			end

			% store the raw data for reasons that have yet to become clear
			rawcal.Vrms(n) = cal.Vrms(n);
			rawcal.unfilt_rms(n, :) = unfilt_rms;
			rawcal.unfilt_dbspl(n, :) = unfilt_dbspl;
			rawcal.filt_rms(n, :) = filt_rms;
			rawcal.filt_dbspl(n, :) = filt_dbspl;

			% save the data for this RMS value
			RMScal(n, 1) = cal.Vrms(n);
			RMScal(n, 2) = mean(unfilt_rms);
			RMScal(n, 3) = mean(unfilt_dbspl);
			RMScal(n, 4) = mean(filt_rms);
			RMScal(n, 5) = mean(filt_dbspl);

			RMScal(n, 6) = std(unfilt_rms);
			RMScal(n, 7) = std(unfilt_dbspl);
			RMScal(n, 8) = std(filt_rms);
			RMScal(n, 9) = std(filt_dbspl);
		end

		% save the rms and level data vectors
		caldata.v_rms = RMScal(:, 1);
		caldata.pa_rms = RMScal(:, 4);
		caldata.dbspl = RMScal(:, 5);
		caldata.pa_rms_std = RMScal(:, 8);
		caldata.dbspl_std = RMScal(:, 9);
		cal.RMScal = RMScal;

		% plots
		subplot(413)
		errorbar(caldata.v_rms, caldata.pa_rms, caldata.pa_rms_std);
		title('Vrms vs. Pa rms');
		xlabel('Output Vrms');
		ylabel('Output Pa rms');
		subplot(414)
		errorbar(caldata.v_rms, caldata.dbspl, caldata.dbspl_std);
		title('Vrms vs. dB SPL');
		xlabel('Output Vrms');
		ylabel('Output dB SPL');
		drawnow

		% compute regression fit to db calibration data
		y = caldata.v_rms;
		X = [ones(size(caldata.pa_rms)) caldata.pa_rms];
		[b, bint, r, rint, stats] = regress(y, X);
		caldata.dbcal.intercept = b(1);
		caldata.dbcal.slope = b(2);
		caldata.dbcal.bint = bint;
		caldata.dbcal.r = r;
		caldata.dbcal.rint = rint;
		caldata.dbcal.Rsq = stats(1);
		caldata.dbcal.F = stats(2);
		caldata.dbcal.p = stats(3);
		caldata.dbcal.err = stats(4);
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Analyze tone data
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% store the individual speaker calibration data
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if SAVE_SPEAKERDATA
			save(speakernetfile, 'caldata');
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% store data for this speaker in the netcal{} array
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		netcal{spkrIndex} = caldata;	
	end

%------------------------------------------------------------------------
% store the calibration data
%------------------------------------------------------------------------
	if SAVE_NETDATA
		save(netcalfile, 'netcal');
	end

%------------------------------------------------------------------------
% clean up
%------------------------------------------------------------------------
	warning on MATLAB:log:logOfZero
