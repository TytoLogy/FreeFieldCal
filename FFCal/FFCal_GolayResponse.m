% FFCal_GolayResponse.m
% 
% Run the Golay Calibration Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created: 2007
% (re)Created: 8 June, 2009
%
% Revisions:
%------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for SPL NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	arraySetTiming(outdev, indev, cal.SweepPeriod, cal.StimDelay,...
									cal.StimDuration, cal.AcqDuration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the start and end bins for the Golay calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	start_bin = ms2bin(cal.StimDelay, indev.Fs);
	if ~start_bin
		start_bin = 1;
	end
	end_bin = ms2bin(cal.StimDelay, indev.Fs) + length(cal.A);
	zerostim = zeros(ms2bin(cal.StimDelay + cal.StimDuration, indev.Fs), 1);
	cal.golaystartbin = start_bin;
	cal.golayendbin = end_bin;
	ADOFFSET = ms2samples(0.1, indev.Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run through the golay stimuli/responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% compute # of points to acquire
	acqpts = ms2samples(cal.AcqDuration, indev.Fs);

	for rep = 1:cal.NCalreps
		disp(sprintf('...Golay (frequency) calibration, rep = %d', rep));

		pause(0.5);
		% play stimulus A of Golay pair;
		ar = arraysingleIO(outdev, indev, zBUS, Astim, outdev.channel, acqpts);
		% play stimulus B of Golay pair;
		pause(0.5);
		br = arraysingleIO(outdev, indev, zBUS, Bstim, outdev.channel, acqpts);

		% store the response data
		speakercal.Ar{rep} = ar(start_bin:end_bin);
		speakercal.Br{rep} = br(start_bin:end_bin);

		% Plot the data
		stim_t = 1000 * (0:Alen-1) ./ outdev.Fs;
		resp_t = 1000 * (0:length(ar)-1) ./ indev.Fs;

		if rep == 1
			% Set up the plots
			figure(stimrespfig1);
			subplot(4, 1, 1)
			stim_a = plot(stim_t, cal.A, 'XDataSource', 'stim_t', 'YDataSource', 'cal.A');
			title('Golay Calibration');
			ylabel('Stim A (V)');
			subplot(4, 1, 2)
			resp_a = plot(resp_t, ar, 'XDataSource', 'resp_t', 'YDataSource', 'ar');
			ylabel('Resp A (V)')
			subplot(4, 1, 3)
			stim_b = plot(stim_t, cal.B, 'XDataSource', 'stim_t', 'YDataSource', 'cal.B');
			ylabel('Stim B (V)');
			subplot(4, 1, 4)
			resp_b = plot(resp_t, br, 'XDataSource', 'resp_t', 'YDataSource', 'br');
			xlabel('Time (msec)');
			ylabel('Resp B (V)');
			drawnow;
			plot_init = 0;
		else
			% refresh plots
			figure(stimrespfig1);
			stim_t = 1000 * (0:Alen-1) ./ outdev.Fs;
			resp_t = 1000 * (0:length(ar)-1) ./ indev.Fs;
			refreshdata(stim_a);
			refreshdata(resp_a);
			refreshdata(stim_b);
			refreshdata(resp_b);
			drawnow;
		end
	end	%reps loop end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, process the Golay Calibration Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('...computing transfer function, speaker [%d, %d]', az, el));
	
	H = cell(cal.NCalreps, 1);
	I = cell(cal.NCalreps, 1);
	Tg = cell(cal.NCalreps, 1);
	Pg = cell(cal.NCalreps, 1);

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
	if cal.NCalreps > 1
		Tgm=mean(Tgm);
		Pgm=mean(Pgm);
	end

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
	% smooth the data - in raw form, it's going to be noisy
	caldata.maginv = smooth(caldata.maginv, SMOOVE_SPAN, SMOOVE_METHOD)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot mean response and trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure(xferfig1);
	c=['r' 'm' 'g' 'k' 'y' 'b'];
	for n=1:cal.NCalreps
		subplot(2, 1, 1)
		hold on 
		plot(Freq, db(Tg{n}), c(n)) 
		hold off
		subplot(2, 1, 2)
		hold on
		plot(Freq, Pg{n}, c(n));
		hold off
	end
	subplot(2, 1, 1)
	hold on
	plot(Freq, db(Tgm), 'k') 
	hold off
	title(sprintf('Transfer Function Az: %d El %d', az, el));
	ylabel('Gain (dB)');
	subplot(2, 1, 2)
	hold on
	plot(Freq, Pgm, 'k');
	hold off
	xlabel('Frequency (Hz)');
	ylabel('Phase');
	set(gcf, 'Name', sprintf('Az: %d El %d Transfer Function', az, el));
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_xferfunction' '.fig'];
	saveas(xferfig1, figfilename)

	figure(xferfig2);
	subplot(3, 1, 1)
	plot(caldata.freq, magdb) 
	ylabel('Mag (db)');
	subplot(3, 1, 2)
	plot(caldata.freq, db(caldata.maginv))
	ylabel('Mag Inverse (db)');
	subplot(3, 1, 3)
	plot(caldata.freq, caldata.phase);
	xlabel('Frequency (Hz)');
	ylabel('Phase');
	drawnow;
	set(gcf, 'Name', sprintf('Az: %d El %d Inverse Function', az, el));
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_inversefunction' '.fig'];
	saveas(xferfig2, figfilename)
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now synthesize a noise and play it to test it out
% select the zero speaker and synthesize a dummy stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% synthesize and window uncalibrated sound
	teststim  = synmononoise_fft(cal.StimDuration, outdev.Fs, Gtest_lo, Gtest_hi, UncalibratedGolayTestRMS, 0);		
	teststim = sin2array(teststim, cal.TestRamp, outdev.Fs);
	
	% synthesize and window calibrated sound
	testcalstim  = synmononoise_fft(cal.StimDuration, outdev.Fs, Gtest_lo, Gtest_hi, CalibratedGolayTestRMS, caldata);		
	testcalstim = sin2array(testcalstim, cal.TestRamp, outdev.Fs);

	% play the test stimuli, get responses
	disp('...playing uncalibrated noise');
	testresp = arraysingleIO(outdev, indev, zBUS, teststim, outdev.channel, ms2samples(cal.AcqDuration, indev.Fs));
	disp('...playing calibrated noise');
	testcalresp = arraysingleIO(outdev, indev, zBUS, testcalstim, outdev.channel, ms2samples(cal.AcqDuration, indev.Fs));

	%filter the test stimuli (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, testresp);
	testcalresp = filter(cal.coeffb, cal.coeffa, testcalresp);
	
	% get the response only
	start_bin = ms2samples(cal.StimDelay, indev.Fs);
	end_bin = start_bin + ms2samples(cal.StimDuration + cal.TestRamp, indev.Fs);
	testresp = testresp(start_bin:end_bin);
	testcalresp = testcalresp(start_bin:end_bin);

	% compute signal levels (not necessary here, but useful diagnostics)
	testresp_dB = dbspl(rms(VtoPa*testresp));
	testcalresp_dB = dbspl(rms(VtoPa*testcalresp));
	
	disp(sprintf('...uncalibrated dB: %.4f (Vrms set at %.4f)',...
							testresp_dB, UncalibratedGolayTestRMS))
	
	disp(sprintf('...calibrated dB: %.4f (Vrms set at %.4f)',...
							testcalresp_dB, CalibratedGolayTestRMS))

	% plot data and response FFT
	fftplot(testresp, indev.Fs, stimrespfig1);
	subplot(311)
	title('Golay Cal: Uncalibrated noise');
	fftplot(testcalresp, indev.Fs, stimrespfig2);
	subplot(311)
	title('Golay Cal: Calibrated noise');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrap it up... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('Golay Frequency/Phase Calibration Finished for speaker [%d, %d]', az, el))

