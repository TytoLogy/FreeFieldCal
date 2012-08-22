warning off MATLAB:log:logOfZero

% this sets up the settings from the data acquisition
run C:\Users\Matlab\Calibration\Data\13-Jan-2008\FFCal_settings_13Jan2008

% specify the input directory
indatadir = 'C:\Users\Matlab\Calibration\Data\13-Jan-2008';

% specify the output directory
outdatadir = [pwd '\NetCal'];
mkdir(outdatadir);

load([indatadir '\ArrayCalData.mat']);
indev = cal.indev;
outdev = cal.outdev;

netcal = cell(Nspeakers, 1);

az = input('azimuth? ');
el = input('elevation? ');

currentspkr = speaker_data('azel', [az el])
az = currentspkr.azimuths
el = currentspkr.elevations

fname = sprintf('%d_%d', az, el);
fname(find(fname=='-')) = 'n';		
speakercalfile = [indatadir '\' fname '.mat']

load(speakercalfile);

caldata.speakerinfo = currentspkr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, process the Golay Calibration Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('Computing Transfer function, speaker az %d el %d ...', az, el));

	for n=1:Nreps
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

	% plot mean response and trials
	c=['r' 'm' 'g' 'k' 'y' 'b'];
	subplot(4, 1, 1),	plot(Freq, db(Tg{n}), c(n));
	subplot(4, 1, 2),	plot(Freq, Pg{n}, c(n));
	for n=1:Nreps
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze the RMS-SPL  Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for n = 1:cal.Nrms
		for r = 1:cal.RMSreps
			splfilt = filter(cal.coeffb, cal.coeffa, speakercal.SPL{n, r});
			unfilt_rms(r) = rms(VtoPa * speakercal.SPL{n, r});
			unfilt_dbspl(r) = dbspl(unfilt_rms(r));
			filt_rms(r) = rms(VtoPa * splfilt);
			filt_dbspl(r) = dbspl(filt_rms(r));
		end

		rawcal.Vrms(n) = cal.Vrms(n);
		rawcal.unfilt_rms(n, :) = unfilt_rms;
		rawcal.unfilt_dbspl(n, :) = unfilt_dbspl;
		rawcal.filt_rms(n, :) = filt_rms;
		rawcal.filt_dbspl(n, :) = filt_dbspl;

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

warning on MATLAB:log:logOfZero


	