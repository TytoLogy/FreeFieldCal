% FFCal_SPLbandwidth.m
%
% Computes the SPL vs. noise bandwidth relationship
%
% Sharad Shanbhag
% sshanbha@aecom.yu.edu

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created: 
%	1 June, 2009
% 
% Revisions:
%  8 June, 2009: rebuilt
%------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the SPLBandwidth Curve Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Running SPL - bandwidth determination');
	pause(0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute # of samples to collect
	in_samples = ms2samples(cal.SPLAcquisition, indev.Fs);
	% start bin for rms signal
	cal.splstart = ms2bin(cal.SPLDelay + cal.SPLRamp, indev.Fs);
	% end bin for rms signal
	cal.splend = cal.splstart + ms2bin(cal.SPLDuration, indev.Fs);
	ndatapts = cal.splend-cal.splstart + 1;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for SPL NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	arraySetTiming(outdev, indev, cal.SPLSweepPeriod, cal.SPLDelay,...
									cal.SPLDuration, cal.SPLAcquisition);
	% Set the sweep count (may not be necessary)
	RPsettag(outdev, 'SwCount', 1);
	RPsettag(indev, 'SwCount', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play calibrated noise at different  bandwidths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	% loop through the different noise bandwidths (at single
	% center frequency) specified in cal.bwfreqs (see FFCal_settings.m)
	disp('SPL-BW calibration')
	
	for rms_index = 1:cal.Nsplbwrms
		
		rmsval = cal.SPL_BWrms(rms_index);
		disp(sprintf('\t RMS = %.4f', rmsval));

		for bw_index = 1:cal.Nbwfreqs
			BW = cal.bwfreqs(bw_index);
			fmin = cal.bwctrfreq - BW / 2;
			fmax = cal.bwctrfreq + BW / 2;
			disp(sprintf('...BW = %.3f', BW));
			
			% preallocate the data array
			speakercal.SPL_BW{rms_index, bw_index} = zeros(cal.SPL_BWreps, ndatapts);
			
			% loop through # of reps
			for rep = 1:cal.SPL_BWreps
				% synthesize and window the calibrated sound &
				% set RMS scale value to cal.SPL_BWrms
				rms_stim = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, rmsval, caldata);
				rms_stim = sin2array(rms_stim, cal.SPLRamp, outdev.Fs);

				% play stimulus, filter the response;
				splresp = arraysingleIO(outdev, indev, zBUS, rms_stim, outdev.channel, in_samples);
				splrespfilt = filter(cal.coeffb, cal.coeffa, splresp(cal.splstart:cal.splend));

				% store the raw response vector
				speakercal.SPL_BW{rms_index, bw_index}(rep, :) = splresp(cal.splstart:cal.splend);

				% plot
				splstim_t = 1000 * (0:length(rms_stim)-1) ./ outdev.Fs;
				splresp_t = 1000 * (0:length(splrespfilt)-1) ./ indev.Fs;

				if bw_index == 1
					% Set up the plots
					figure(stimrespfig1);
					subplot(2, 1, 1)
					splstimplot = plot(splstim_t, rms_stim, 'XDataSource', 'splstim_t', 'YDataSource', 'rms_stim');
					title('SPL Calibration');	ylabel('SPL Stim (V)');
					subplot(2, 1, 2)
					splrespplot = plot(splresp_t, splrespfilt, 'XDataSource', 'splresp_t', 'YDataSource', 'splrespfilt');
					xlabel('Time (msec)');		ylabel('Resp (V)')
					bwplot_init = 0;
				else
					figure(stimrespfig1);
					% refresh plots
					refreshdata(splstimplot);
					refreshdata(splrespplot);
					drawnow;
				end

			end	% end of REP loop
		end	% end of BW loop
	end % end of RMSindex loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze and store the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tmp = zeros(cal.Nsplbwrms, cal.Nbwfreqs);
	caldata.bw.rms_mean = tmp;
	caldata.bw.rms_std = tmp;
	caldata.bw.db_mean = tmp;
	caldata.bw.db_std = tmp;

	% loop through the # of Vrms values for the BW/SPL data (cal.Nsplbwrms)
	for rms_index = 1:cal.Nsplbwrms
		
		% preallocate some arrays
		tmp = zeros(cal.SPL_BWreps, 1);
		filt_rms = tmp;
		filt_dbspl = tmp;
		
		% loop through the # of bandwidths tested (cal.Nbwfreqs)
		for bw_index = 1:cal.Nbwfreqs
			% get the data for this bandwidth
			tmpdata = speakercal.SPL_BW{rms_index, bw_index};
			
			% loop through reps (cal.SPL_BWreps)
			for r = 1:cal.SPL_BWreps
				filt_rms(r) = rms( VtoPa * filter(cal.coeffb, cal.coeffa, tmpdata(r, :)) );
				filt_dbspl(r) = dbspl(filt_rms(r));
			end
			
			% compute mean and std deviation for data (if more than 1 rep)
			if cal.SPL_BWreps > 1
				caldata.bw.rms_mean(rms_index, bw_index) = mean(filt_rms);
				caldata.bw.rms_std(rms_index, bw_index) = std(filt_rms);
				caldata.bw.db_mean(rms_index, bw_index) = mean(filt_dbspl);
				caldata.bw.db_std(rms_index, bw_index) = std(filt_dbspl);
			else
				caldata.bw.rms_mean(rms_index, bw_index) = filt_rms;
				caldata.bw.db_mean(rms_index, bw_index) = filt_dbspl;
			end
		end		% end of bw_index loop
	end		% end of rms_index loop

	caldata.bw.bwfreqs = cal.bwfreqs;
	% compute max_bw difference
	caldata.bw.maxdb = max(caldata.bw.db_mean(cal.Nsplbwrms, cal.Nbwfreqs));
	% compute the correction (in db) vs. bandwidth
	caldata.bw.correction = caldata.bw.maxdb - caldata.bw.db_mean(cal.Nsplbwrms, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       BW_DATAFIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create plot
	figure(bw_datafig);
	% create the legend
	syms = '.ox+*sdv^<>ph';
	for n = 1:cal.Nsplbwrms
		legstr{n} = num2str(cal.SPL_BWrms(n));
	end

	subplot(121)
	plot(caldata.bw.bwfreqs, caldata.bw.rms_mean(1, :), [syms(1) '-'])
	if cal.Nsplbwrms > 1
		hold on
		for f = 2:cal.Nsplbwrms
			plot(caldata.bw.bwfreqs, caldata.bw.rms_mean(f, :), [syms(f) '-'])
		end
		hold off
	end
	titlestr{1} = sprintf('RMS(Pa) vs. BW (Hz), az: %d el: %d', az, el);
	titlestr{2} = sprintf('Fctr= %d, Vrms = %.2f', ...
							cal.bwctrfreq, cal.SPL_BWrms(1));
	title(titlestr)
	xlabel('Output Signal Bandwidth (Hz)');
	ylabel('Output Pa rms');
	grid on
	legend(legstr);
	
	subplot(122)
	plot(caldata.bw.bwfreqs, caldata.bw.db_mean(1, :), [syms(1) '-'])
	if cal.Nsplbwrms > 1
		hold on
		for f = 2:cal.Nsplbwrms
			plot(caldata.bw.bwfreqs, caldata.bw.db_mean(f, :), [syms(f) '-'])
		end
		hold off
	end
	title('dB SPL vs. Output Signal Bandwidth (Hz) ');
	xlabel('Output Signal Bandwidth (Hz)');
	ylabel('Output dB SPL');
	grid on
	legend(legstr);

	set(gcf, 'Name', 'Pa vs. Bandwidth');
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_Pa_vs_BW' '.fig'];
	saveas(bw_datafig, figfilename)

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       BW_SEMILOG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	figure(bw_semilog);
	subplot(121)
	semilogx(caldata.bw.bwfreqs, caldata.bw.rms_mean(1, :), [syms(1) '-'])
	if cal.Nsplbwrms > 1
		hold on
		for f = 2:cal.Nsplbwrms
			semilogx(caldata.bw.bwfreqs, caldata.bw.rms_mean(f, :), [syms(f) '-'])
		end
		hold off
	end
	titlestr{1} = sprintf('RMS(Pa) vs. BW (Hz), az: %d el: %d', az, el);
	titlestr{2} = sprintf('Fctr= %d, Vrms = %.2f', ...
							cal.bwctrfreq, cal.SPL_BWrms(1));
	title(titlestr)
	xlabel('Output Signal Bandwidth (log_1_0(Hz))');
	ylabel('Output Pa rms');
	grid on
	legend(legstr);
	
	subplot(122)
	semilogx(caldata.bw.bwfreqs, caldata.bw.db_mean(1, :), [syms(1) '-'])
	if cal.Nsplbwrms > 1
		hold on
		for f = 2:cal.Nsplbwrms
			semilogx(caldata.bw.bwfreqs, caldata.bw.db_mean(f, :), [syms(f) '-'])
		end
		hold off
	end
	title('dB SPL vs. Output Signal Bandwidth (Hz) ');
	xlabel('Output Signal Bandwidth (log_1_0(Hz))');
	ylabel('Output dB SPL');
	grid on
	legend(legstr);
	set(gcf, 'Name', 'dB SPL vs. Bandwidth');
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_dB_vs_BW' '.fig'];
	saveas(bw_semilog, figfilename)
	
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now synthesize and play uncorrected 
% broadband and narrowband noise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Broadband and Narrowband noise, no BW correction...')

	% get the stimulus intensity level scaling factor from the 
	% calibration data
	stimlevel_pa = dbspl2pa(cal.BWtestlevel);
	scale = get_scale(stimlevel_pa, caldata.v_rms, caldata.pa_rms);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%         BROADBAND
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% synthesize the broadband stimulus
	fmin = cal.BWbroad.fctr(1) - cal.BWbroad.bw(1) / 2;
	fmax = cal.BWbroad.fctr(1) + cal.BWbroad.bw(1) / 2;
	bbstim  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
	% ramp the stimulus on and off for smoothness' sake...
	bbstim = sin2array(bbstim, cal.TestRamp, outdev.Fs);

	% play the test broadband stimulus
	bbresp = arraysingleIO(outdev, indev, zBUS, bbstim, outdev.channel, in_samples);

	%filter the test stimulus (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, bbresp);
	testresp = testresp(cal.splstart:cal.splend);
	testdbspl = dbspl(rms(VtoPa*testresp));

	disp(sprintf('\tDesired bb stim dB SPL: %f', cal.BWtestlevel));
	disp(sprintf('\tMeasured test stim dB SPL: %f', testdbspl));

	fftplot(testresp, indev.Fs, stimrespfig1);
	subplot(311)
	title('Freq Resp, Uncorrected BB RMS Test Stimulus');
	set(gcf, 'Name', 'Uncorrected FR, BB Stimulus');
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_BBuncorr_rmstest' '.fig'];
	saveas(stimrespfig1, figfilename)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%         NARROWBAND
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% synthesize the narrow band stimulus
	fmin = cal.BWnarrow.fctr(1) - cal.BWnarrow.bw(1) / 2;
	fmax = cal.BWnarrow.fctr(1) + cal.BWnarrow.bw(1) / 2;
	nbstim  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
	% ramp the stimulus on and off for smoothness' sake...
	nbstim = sin2array(nbstim, cal.TestRamp, outdev.Fs);

	% play the test broadband stimulus
	nbresp = arraysingleIO(outdev, indev, zBUS, nbstim, outdev.channel, in_samples);

	%filter the test stimulus (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, nbresp);
	testresp = testresp(cal.splstart:cal.splend);
	testdbspl = dbspl(rms(VtoPa*testresp));

	disp(sprintf('\tDesired nb stim dB SPL: %f', cal.BWtestlevel));
	disp(sprintf('\tMeasured nb test stim dB SPL: %f', testdbspl));

	fftplot(testresp, indev.Fs, stimrespfig2);
	subplot(311)
	title('Freq Resp, Uncorrected NB RMS Test Stimulus');
	set(gcf, 'Name', 'Uncorrected FR, NB Stimulus');
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_NBuncorr_rmstest' '.fig'];
	saveas(stimrespfig1, figfilename)
	
	pause(1);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test scaling  using get_bw_scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Testing BW correction SPL factor...')
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%         BROADBAND
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% synthesize the broadband stimulus
	fmin = cal.BWbroad.fctr(1) - cal.BWbroad.bw(1) / 2;
	fmax = cal.BWbroad.fctr(1) + cal.BWbroad.bw(1) / 2;
	% get the stimulus intensity level scaling factor from the 
	% calibration data
	scale = get_bw_scale_nowarning(cal.BWtestlevel, caldata, [fmin fmax]);

	[bbstim, bbmag, bbphase]  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
	% ramp the stimulus on and off for smoothness' sake...
	bbstim = sin2array(bbstim, cal.TestRamp, outdev.Fs);

	% play the test broadband stimulus
	bbresp = arraysingleIO(outdev, indev, zBUS, bbstim, outdev.channel, in_samples);

	%filter the test stimulus (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, bbresp);
	testresp = testresp(cal.splstart:cal.splend);
	testdbspl = dbspl(rms(VtoPa*testresp));

	disp(sprintf('\tDesired bb stim dB SPL: %f', cal.BWtestlevel));
	disp(sprintf('\tMeasured test stim dB SPL: %f', testdbspl));

	fftplot(testresp, indev.Fs, stimrespfig1);
	subplot(311)
	title('Frequency Response, BB RMS Test Stimulus');
	set(gcf, 'Name', 'Freq Resp, BB Stimulus');
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_BBrmstest' '.fig'];
	saveas(stimrespfig1, figfilename)
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%         NARROWBAND
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% synthesize the narrow band stimulus
	fmin = cal.BWnarrow.fctr(1) - cal.BWnarrow.bw(1) / 2;
	fmax = cal.BWnarrow.fctr(1) + cal.BWnarrow.bw(1) / 2;
	% get the stimulus intensity level scaling factor from the 
	% calibration data
	scale = get_bw_scale_nowarning(cal.BWtestlevel, caldata, [fmin fmax]);
	[nbstim, nbmag, nbphase]  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
	% ramp the stimulus on and off for smoothness' sake...
	nbstim = sin2array(nbstim, cal.TestRamp, outdev.Fs);

	% play the test narrowband stimulus
	nbresp = arraysingleIO(outdev, indev, zBUS, nbstim, outdev.channel, in_samples);

	%filter the test stimulus (eliminate the nasty low frequency crap)
	testresp = filter(cal.coeffb, cal.coeffa, nbresp);
	testresp = testresp(cal.splstart:cal.splend);
	testdbspl = dbspl(rms(VtoPa*testresp));

	disp(sprintf('\tDesired nb stim dB SPL: %f', cal.BWtestlevel));
	disp(sprintf('\tMeasured nb test stim dB SPL: %f', testdbspl));

	fftplot(testresp, indev.Fs, stimrespfig2);
	subplot(311)
	title('Frequency Response, NB RMS Test Stimulus');
	set(gcf, 'Name', 'Freq Resp, NB Stimulus');
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_NBrmstest' '.fig'];
	saveas(stimrespfig2, figfilename)
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrap it up... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('BW - SPL Calibration Finished for speaker [%d, %d]', az, el))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FFCal_SPLbandwidth.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%