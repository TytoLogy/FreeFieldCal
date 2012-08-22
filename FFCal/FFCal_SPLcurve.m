%------------------------------------------------------------------------
% FFCal_SPLcurve.m
%------------------------------------------------------------------------
% 
% Determines Vrms -> Pa (& dB SPL) conversion for broadband noise
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sharad.shanbhag@einstein.yu.edu
%------------------------------------------------------------------------
% Created: 30 Oct 2007
%
% Revisions:
%	10 Jan 2008: added plotting feature
%	1 June, 2009: cleaned up a bit, added SPL-specific timing/duration code
%  8 June, 2009: rebuilt
%	18 March, 2010 (SJS): updated comments
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the SPL Curve Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Running SPL curve determination');
	pause(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute # of samples to collect
	in_samples = ms2samples(cal.SPLAcquisition, indev.Fs);
	% start bin for rms signal
	cal.splstart = ms2bin(cal.SPLDelay + cal.SPLRamp, indev.Fs);
	% end bin for rms signal
	cal.splend = cal.splstart + ms2bin(cal.SPLDuration, indev.Fs);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for SPL NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	arraySetTiming(outdev, indev, cal.SPLSweepPeriod, cal.SPLDelay,...
									cal.SPLDuration, cal.SPLAcquisition);
	% Set the sweep count (may not be necessary)
	RPsettag(outdev, 'SwCount', 1);
	RPsettag(indev, 'SwCount', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play calibrated noise to compute SPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for rms_index = 1:cal.Nrms
		disp(sprintf('...intensity calibration, RMS = %.3f', cal.Vrms(rms_index)));

		for rep = 1:cal.RMSreps
			% synthesize calibrated rms stimulus
			rms_stim = synmononoise_fft(cal.SPLDuration, outdev.Fs, cal.splminfreq, cal.splmaxfreq, cal.Vrms(rms_index), caldata);
			rms_stim = sin2array(rms_stim, cal.SPLRamp, outdev.Fs);

			% play stimulus;
			splresp = arraysingleIO(outdev, indev, zBUS, rms_stim, outdev.channel, in_samples);
			splrespfilt = filter(cal.coeffb, cal.coeffa, splresp(cal.splstart:cal.splend));
			speakercal.SPL{rms_index, rep} = splresp(cal.splstart:cal.splend);

			tmp = filter(cal.coeffb, cal.coeffa, speakercal.SPL{rms_index, rep});
			tmpdb = dbspl(rms(VtoPa * tmp));
			disp(sprintf('\t...rep %d \t %.2f dB SPL', rep, tmpdb));

			if rms_index == 1
				% Set up the plots
				figure(stimrespfig1);
				splstim_t = 1000 * (0:length(rms_stim)-1) ./ outdev.Fs;
				splresp_t = 1000 * (0:length(splrespfilt)-1) ./ indev.Fs;			
				subplot(2, 1, 1)
				splstimplot = plot(splstim_t, rms_stim, 'XDataSource', 'splstim_t', 'YDataSource', 'rms_stim');
				title('SPL Calibration');	ylabel('SPL Stim (V)');
				subplot(2, 1, 2)
				splrespplot = plot(splresp_t, splrespfilt, 'XDataSource', 'splresp_t', 'YDataSource', 'splrespfilt');
				xlabel('Time (msec)');
				ylabel('Resp (V)')
				set(gcf, 'Name', 'SPL Curve Stim & Response')
				drawnow;
				
			else
				% refresh plots
				figure(stimrespfig1);
				refreshdata(splstimplot);
				refreshdata(splrespplot);
				drawnow;
			end
		end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze the Data and store the calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% pre-allocate storage for RMS calibration data
	RMScal = zeros(cal.Nrms, 9);
	
	tmp = zeros(cal.RMSreps, 1);
	for n = 1:cal.Nrms
		unfilt_rms = tmp;
		unfilt_dbspl = tmp;
		filt_rms = tmp;
		filt_dbspl = tmp;

		for r = 1:cal.RMSreps
			splfilt = filter(cal.coeffb, cal.coeffa, speakercal.SPL{n, r});
			unfilt_rms(r) = rms(VtoPa * speakercal.SPL{n, r});
			unfilt_dbspl(r) = dbspl(unfilt_rms(r));
			filt_rms(r) = rms(VtoPa * splfilt);
			filt_dbspl(r) = dbspl(filt_rms(r));
		end

		RMScal(n, 1) = cal.Vrms(n);
		if cal.RMSreps > 1
			RMScal(n, 2) = mean(unfilt_rms);
			RMScal(n, 3) = mean(unfilt_dbspl);
			RMScal(n, 4) = mean(filt_rms);
			RMScal(n, 5) = mean(filt_dbspl);
			RMScal(n, 6) = std(unfilt_rms);
			RMScal(n, 7) = std(unfilt_dbspl);
			RMScal(n, 8) = std(filt_rms);
			RMScal(n, 9) = std(filt_dbspl);
		else
			RMScal(n, 2) = (unfilt_rms);
			RMScal(n, 3) = (unfilt_dbspl);
			RMScal(n, 4) = (filt_rms);
			RMScal(n, 5) = mean(filt_dbspl);
		end
	end

	speakercal.RMScal = RMScal;

	caldata.v_rms = RMScal(:, 1);
	caldata.pa_rms = RMScal(:, 4);
	caldata.dbspl = RMScal(:, 5);
	caldata.pa_rms_std = RMScal(:, 8);
	caldata.dbspl_std = RMScal(:, 9);
	
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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot the response curve
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure(splcurvefig);
	subplot(121)
	errorbar(caldata.v_rms, caldata.pa_rms, caldata.pa_rms_std, '.-');
	title(sprintf('Vrms vs. Pa rms, az: %d el: %d', az, el));
	xlabel('Output Vrms');
	ylabel('Output Pa rms');
	subplot(122)
	errorbar(caldata.v_rms, caldata.dbspl, caldata.dbspl_std, '.-');
	title('Vrms vs. dB SPL');
	xlabel('Output Vrms');
	ylabel('Output dB SPL');
	set(gcf, 'Name', 'Output Level vs. Signal Vrms')
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_SPLvsVrms' '.fig'];
	saveas(splcurvefig, figfilename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now synthesize a noise and play it to test it out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% get the stimulus intensity level scaling factor from the 
	% calibration data
	stimlevel_pa = dbspl2pa(cal.SPLtestlevel_db);
	scale = get_scale(stimlevel_pa, caldata.v_rms, caldata.pa_rms);

	[teststim, tmag, tphase]  = synmononoise_fft(cal.SPLDuration, outdev.Fs, cal.splminfreq, cal.splmaxfreq, scale, caldata);		
	% ramp the stimulus on and off for smoothness' sake...
	S = sin2array(teststim, cal.TestRamp, outdev.Fs);

	% play the test stimulus
	testresp = arraysingleIO(outdev, indev, zBUS, teststim, outdev.channel, ms2samples(cal.AcqDuration, indev.Fs));

	%filter the test stimulus (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, testresp);
	t1 = ms2samples(cal.StimDelay, indev.Fs);
	t2 = t1 + ms2samples(cal.SPLDuration + cal.TestRamp, indev.Fs);
	testresp = testresp(t1:t2);
	testdbspl = dbspl(rms(VtoPa*testresp));

	disp(sprintf('\tDesired test stim dB SPL: %.2f', cal.SPLtestlevel_db));
	disp(sprintf('\tMeasured test stim dB SPL: %.2f', testdbspl));
	disp(sprintf('\t\t\t error: %.4f dB SPL', testdbspl - cal.SPLtestlevel_db));

	fftplot(testresp, indev.Fs, stimrespfig2);
	subplot(311)
	title('Frequency Response, RMS Test Stimulus');
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrap it up... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('RMS - SPL Calibration Finished for speaker [%d, %d]', az, el))



