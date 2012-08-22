% tester.m
%
% Gets calibration from B&K mic
%
% Sharad Shanbhag
% sshanbha@aecom.yu.edu

%
%	2 June, 2009: created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	FFCal_settings;
	FFCal_filesettings;

	ADOFFSET = ms2samples(0.1, indev.Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the TDT circuits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	FFCal_tdtinit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get a highpass filter for processing the RMS data
%  (needed to wait until Fs is known for indev)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Nyquist frequency
	cal.fnyq = indev.Fs/2;
	
	[cal.coeffb, cal.coeffa] = butter(cal.forder, cal.fcutoff/cal.fnyq, 'high');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute # of samples to collect
	in_samples = ms2samples(cal.SPLAcquisition, indev.Fs);

	% start bin for rms signal
	cal.splstart = ms2bin(cal.SPLDelay + cal.SPLRamp, indev.Fs);
	% end bin for rms signal
	cal.splend = cal.splstart + ms2bin(cal.SPLDuration, indev.Fs);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for SPL NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Set the total sweep period time
	RPsettag(outdev, 'SwPeriod', ms2samples(cal.SPLSweepPeriod, outdev.Fs));
	RPsettag(indev, 'SwPeriod', ms2samples(cal.SPLSweepPeriod, indev.Fs));
	% Set the sweep count (may not be necessary)
	RPsettag(outdev, 'SwCount', 1);
	RPsettag(indev, 'SwCount', 1);
	% Set the Stimulus Delay
	RPsettag(outdev, 'StimDelay', ms2samples(cal.SPLDelay, outdev.Fs));
	% Set the Stimulus Duration
	RPsettag(outdev, 'StimDur', ms2samples(cal.SPLDuration, outdev.Fs));
	% Set the length of time to acquire data
	RPsettag(indev, 'AcqDur', ms2samples(cal.SPLAcquisition, indev.Fs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play calibrated noise at different  bandwidths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% loop through the different noise bandwidths (at single
	% center frequency) specified in cal.bwfreqs (see FFCal_settings.m)
	for bw_index = 1:length(cal.bwfreqs)
		BW = cal.bwfreqs(bw_index);
		
		disp(sprintf('...SPL-BW calibration, BW = %.3f', BW));
		
		fmin = cal.bwctrfreq - BW / 2;
		fmax = cal.bwctrfreq + BW / 2;

		% loop through # of reps
		for rep = 1:cal.SPL_BWreps
			
			% synthesize and window the calibrated sound &
			% set RMS scale value to cal.SPL_BWrms
			rms_stim = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, cal.SPL_BWrms, caldata);
			rms_stim = sin2array(rms_stim, cal.SPLRamp, outdev.Fs);

			% play stimulus, filter the response;
			[splresp, index] = arraysingleIO(outdev, indev, zBUS, rms_stim, outdev.channel, in_samples);
			splrespfilt = filter(cal.coeffb, cal.coeffa, splresp(cal.splstart:cal.splend));
			% store the response vector
			speakercal.SPL_BW{bw_index, rep} = splresp(cal.splstart:cal.splend);

			% plot
			splstim_t = 1000 * [0:length(rms_stim)-1] ./ outdev.Fs;
			splresp_t = 1000 * [0:length(splrespfilt)-1] ./ indev.Fs;

			if bwplot_init
				% Set up the plots
				splfigure = figure;
				subplot(2, 1, 1)
				splstimplot = plot(splstim_t, rms_stim, 'XDataSource', 'splstim_t', 'YDataSource', 'rms_stim');
				title('SPL Calibration');	ylabel('SPL Stim (V)');
				subplot(2, 1, 2)
				splrespplot = plot(splresp_t, splrespfilt, 'XDataSource', 'splresp_t', 'YDataSource', 'splrespfilt');
				xlabel('Time (msec)');		ylabel('Resp (V)')
				bwplot_init = 0;
			else
				figure(splfigure);
				% refresh plots
				refreshdata(splstimplot);
				refreshdata(splrespplot);
				drawnow;
			end
			
		end	% end of REP loop
	end	% end of BW loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze and store the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for n = 1:length(cal.bwfreqs)
		for r = 1:cal.SPL_BWreps
			splfilt = filter(cal.coeffb, cal.coeffa, speakercal.SPL_BW{n, r});
			unfilt_rms(r) = rms(VtoPa * speakercal.SPL_BW{n, r});
			unfilt_dbspl(r) = dbspl(unfilt_rms(r));
			filt_rms(r) = rms(VtoPa * splfilt);
			filt_dbspl(r) = dbspl(filt_rms(r));
		end

		% compute averages
		RMSBWcal(n, 1) = cal.bwfreqs(n);
		RMSBWcal(n, 2) = mean(unfilt_rms);
		RMSBWcal(n, 3) = mean(unfilt_dbspl);
		RMSBWcal(n, 4) = mean(filt_rms);
		RMSBWcal(n, 5) = mean(filt_dbspl);

		RMSBWcal(n, 6) = std(unfilt_rms);
		RMSBWcal(n, 7) = std(unfilt_dbspl);
		RMSBWcal(n, 8) = std(filt_rms);
		RMSBWcal(n, 9) = std(filt_dbspl);
	end

	caldata.bw.bwfreqs = RMSBWcal(:, 1);
	caldata.bw.pa_rms = RMSBWcal(:, 4);
	caldata.bw.dbspl = RMSBWcal(:, 5);
	caldata.bw.pa_rms_std = RMSBWcal(:, 8);
	caldata.bw.dbspl_std = RMSBWcal(:, 9);

	cal.RMSBWcal = RMSBWcal;

	bw_datafig = figure;
	subplot(121)
	errorbar(caldata.bw.bwfreqs, caldata.bw.pa_rms, caldata.bw.pa_rms_std, '.-');
	title('bw vs. Pa rms');
	xlabel('Output Signal Bandwidth (Hz)');
	ylabel('Output Pa rms');
	subplot(122)
	errorbar(caldata.bw.bwfreqs, caldata.bw.dbspl, caldata.bw.dbspl_std, '.-');
	title('Output Signal Bandwidth (Hz) vs. dB SPL');
	xlabel('Output Signal Bandwidth (Hz)');
	ylabel('Output dB SPL');
	
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now synthesize a broadband and narrowband noise 
% and play it to test it out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% get the stimulus intensity level scaling factor from the 
	% calibration data
	stimlevel_pa = dbspl2pa(cal.BWtestlevel);
	scale = get_scale(stimlevel_pa, caldata.v_rms, caldata.pa_rms);

	% synthesize the broadband stimulus
	fmin = cal.BWbroad.fctr(1) - cal.BWbroad.bw(1) / 2;
	fmax = cal.BWbroad.fctr(1) + cal.BWbroad.bw(1) / 2;
	[bbstim, bbmag, bbphase]  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
	% ramp the stimulus on and off for smoothness' sake...
	bbstim = sin2array(bbstim, cal.TestRamp, outdev.Fs);

	% play the test broadband stimulus
	[bbresp, index] = arraysingleIO(outdev, indev, zBUS, bbstim, outdev.channel, in_samples);

	%filter the test stimulus (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, bbresp);
	testresp = testresp(t1:t2);
	testdbspl = dbspl(rms(VtoPa*testresp));

	disp(sprintf('\tDesired bb stim dB SPL: %f', cal.BWtestlevel));
	disp(sprintf('\tMeasured test stim dB SPL: %f', testdbspl));

	bbplotfig = figure;
	fftplot(testresp, indev.Fs, bbplotfig);
	subplot(311)
	title('Frequency Response, BB RMS Test Stimulus');

	% synthesize the narrow band stimulus
	fmin = cal.BWnarrow.fctr(1) - cal.BWnarrow.bw(1) / 2;
	fmax = cal.BWnarrow.fctr(1) + cal.BWnarrow.bw(1) / 2;
	[nbstim, nbmag, nbphase]  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
	% ramp the stimulus on and off for smoothness' sake...
	nbstim = sin2array(nbstim, cal.TestRamp, outdev.Fs);

	% play the test broadband stimulus
	[nbresp, index] = arraysingleIO(outdev, indev, zBUS, nbstim, outdev.channel, in_samples);

	%filter the test stimulus (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, nbresp);
	testresp = testresp(t1:t2);
	testdbspl = dbspl(rms(VtoPa*testresp));

	disp(sprintf('\tDesired nb stim dB SPL: %f', cal.BWtestlevel));
	disp(sprintf('\tMeasured nb test stim dB SPL: %f', testdbspl));

	nbplotfig = figure;
	fftplot(testresp, indev.Fs, nbplotfig);
	subplot(311)
	title('Frequency Response, NB RMS Test Stimulus');
	
	disp('Press key to continue...')
	pause

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play narrowband noise at different Fctr (if specified) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	nBW = length(cal.BWnarrow.fctr);
	if nBW > 1

		% get the stimulus intensity level scaling factor from the 
		% calibration data
		stimlevel_pa = dbspl2pa(cal.BWtestlevel);
		scale = get_scale(stimlevel_pa, caldata.v_rms, caldata.pa_rms);
		
		BWctrspl = zeros(nBW, cal.SPL_BWreps);
		
		disp(['bw scale factor = ' num2str(scale)]);

		ctrplotfig = figure;

		for n = 1:nBW
			fctr = cal.BWnarrow.fctr(n);
			bw = cal.BWnarrow.bw(1);
			
			% synthesize the broadband stimulus
			fmin = fctr - bw / 2;
			fmax = fctr + bw / 2;
			[stim, mag, phase]  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
			% ramp the stimulus on and off for smoothness' sake...
			stim = sin2array(stim, cal.TestRamp, outdev.Fs);

			for r = 1:cal.SPL_BWreps
				% play the test narrowband stimulus
				[resp, index] = arraysingleIO(outdev, indev, zBUS, stim, outdev.channel, in_samples);

				%filter the test stimulus (eliminate the nasty low frequency crap)
				testresp = filter(cal.coeffb, cal.coeffa, resp);
				testresp = testresp(t1:t2);
				testdbspl = dbspl(rms(VtoPa*testresp));

				disp(sprintf('\t rep %d \t Desired dB SPL: %f \t Measured dB SPL: %f',...
									r, cal.BWtestlevel, testdbspl));
				
				BWctrspl(n, r) = testdbspl;

				fftplot(testresp, indev.Fs, ctrplotfig);
				subplot(311)
				title('Frequency Response, NB CTR');
			end
		end

		% compute mean
		caldata.BWctrfreqs = cal.BWnarrow.fctr;
		caldata.BWctrdB = mean(BWctrspl');
		caldata.BWctrdB_std = std(BWctrspl');
		
		% plot
		figure
		errorbar(caldata.BWctrfreqs, caldata.BWctrdB, caldata.BWctrdB_std, 'b.-')
		xlabel('Fctr (Hz)')
		ylabel('dB SPL');
		title('dB vs. Narrow Band Fctr')

	end % END of variable ctr freq IF
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrap it up... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('BW - SPL Calibration Finished for speaker [%d, %d]', az, el))
	disp('Press enter to continue...')
	pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exit gracefully (remove serial objects, etc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	FFCal_exit;
	
