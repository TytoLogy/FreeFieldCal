function netcal = process_ffcal_data(varargin)
%netcal = process_ffcal_data(varargin)
%
% processes free-field speaker calibration data collected using FFCal
% utility.  Output is in the form of the netcal{} array that is used by 
% programs to calibrate the output signals played through the free-field
% speaker array.
% 
% Input Arguments:
% 
% Output Arguments:
% 	netcal		net calibration data, stored in cell array of structures
%
% See also: FFCal.m, get_bw_scale
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created: 10 June, 2009, from ProcessNetcal.m script
% 
% Revisions:
%	18 June, 2009 (SJS): some more work on reading in FFCal_settings.m
%								from the raw data directory
%------------------------------------------------------------------------

% turn of log of Zero warning message
warning off MATLAB:log:logOfZero

SAVE_SPEAKERDATA = 0;
SAVE_NETDATA = 1;

% set up smoothing for the transfer function curves (raw
% golay-code extracted curves are very very noisy/jagged)
SMOOVE = 1;
SMOOVE_METHOD = 'rloess';
SMOOVE_SPAN = 0.1;

if ~nargin
	% specify the input directory and input data name
	indatadir = 'W:\Sharad\Calibration\FreeFieldRig\Array\Raw';
	indataname = '09-Jun-2009';

	% specify the output directory
	outdatapath = 'W:\Sharad\Calibration\FreeFieldRig\Array\Processed';
end

indatapath = [indatadir filesep indataname];
outdatadir = [outdatapath filesep indataname filesep 'NetCal-' indataname];

currentdir = pwd;
cd(indatapath);
run('FFCal_settings');
cd(currentdir);

if ~exist(outdatadir, 'dir')
	mkdir(outdatadir);
end

netcalfile = [outdatadir filesep 'netcal-' indataname '.mat'];

load([indatapath filesep 'ArrayCalData.mat']);

indev = cal.indev;
outdev = cal.outdev;

% pre-allocate some arrays
netcal = cell(cal.Nspeakers, 1);
H = cell(cal.NCalreps, 1);
I = H;
Tg = H;
Pg = H;

unfilt_rms = zeros(cal.RMSreps, 1);
unfilt_dbspl = zeros(cal.RMSreps, 1);
filt_rms = zeros(cal.RMSreps, 1);
filt_dbspl = zeros(cal.RMSreps, 1);

for spkrIndex = 1:cal.Nspeakers
	currentspkr = speaker_data('speakernum', spkrIndex)
	az = currentspkr.azimuths
	el = currentspkr.elevations

	fname = sprintf('%d_%d', az, el);
	fname(find(fname=='-')) = 'n';		
	speakercalfile = [indatadir '\' indataname '\' fname '.mat']
	speakernetfile = [outdatadir '\netcal-' fname '.mat'];

	load(speakercalfile);
	
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
	% Analyze the RMS-SPL  Data 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		RMScal = zeros(cal.Nrms, 9);
		
		for n = 1:cal.Nrms
			for r = 1:cal.RMSreps
				splfilt = filter(cal.coeffb, cal.coeffa, speakercal.SPL{n, r});
				unfilt_rms(r) = rms(VtoPa * speakercal.SPL{n, r});
				unfilt_dbspl(r) = dbspl(unfilt_rms(r));
				filt_rms(r) = rms(VtoPa * splfilt);
				filt_dbspl(r) = dbspl(filt_rms(r));
			end

% 			rawcal.Vrms(n) = cal.Vrms(n);
% 			rawcal.unfilt_rms(n, :) = unfilt_rms;
% 			rawcal.unfilt_dbspl(n, :) = unfilt_dbspl;
% 			rawcal.filt_rms(n, :) = filt_rms;
% 			rawcal.filt_dbspl(n, :) = filt_dbspl;

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

		caldata.v_rms = RMScal(:, 1);
		caldata.pa_rms = RMScal(:, 4);
		caldata.dbspl = RMScal(:, 5);
		caldata.pa_rms_std = RMScal(:, 8);
		caldata.dbspl_std = RMScal(:, 9);

		cal.RMScal = RMScal;

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
	% store the calibration data
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if SAVE_SPEAKERDATA
			save(speakernetfile, 'caldata');
		end
		netcal{spkrIndex} = caldata;	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store the calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SAVE_NETDATA
	save(netcalfile, 'netcal');
end

warning on MATLAB:log:logOfZero

