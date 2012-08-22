%------------------------------------------------------------------------
% FFCal_ToneSPLResponse.m
%------------------------------------------------------------------------
%
% Using calibration from B&K mic
% Play tones to check for proper calibration
%
%------------------------------------------------------------------------
% Part of FFCal program
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created:
%	8 June, 2009
%
% Revisions:
%	28 October, 2009 (SJS): more documentation
%------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play tones to check for proper calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Running Tone SPL Response...');
	pause(0.5);

	% compute # of samples to collect
	in_samples = ms2samples(cal.ToneSPLAcquisition, indev.Fs);
	% start bin for tone signal
	tonestart = ms2bin(cal.ToneSPLDelay + cal.ToneSPLRamp, indev.Fs);
	% end bin for tone signal
	toneend = tonestart + ms2bin(cal.ToneSPLDuration, indev.Fs);
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for TONES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	status = arraySetTiming(outdev, indev, cal.ToneSPLSweepPeriod, cal.ToneSPLDelay,...
									cal.ToneSPLDuration, cal.ToneSPLAcquisition);
	% Set the sweep count (may not be necessary)
	RPsettag(outdev, 'SwCount', 1);
	RPsettag(indev, 'SwCount', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play calibrated tones to compute SPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% open a figure for plotting
	figure(stimrespfig1);
	
	% loop through the # of rms values to test
	for rms_index = 1:cal.Ntonesplrms
		% get the rms value to test
		tonerms = cal.ToneSPLRMS(rms_index);
		
		disp(sprintf('...tone intensity calibration, RMS = %.3f', tonerms));

		% loop through the # of tone frequencies to test
		for tone_index = 1:cal.Ntonesplfreqs
			% get the freq to test
			tonefreq = cal.ToneSPLFreqs(tone_index);
			disp(sprintf('\t Freq = %d', tonefreq));

			% generate a *frequency calibrated* tone with the specified 
			% frequency and RMS value and then apply sin^2 window
			tonestim = synmonosine(cal.ToneSPLDuration, outdev.Fs, tonefreq, tonerms, caldata);
			tonestim = sin2array(tonestim, cal.ToneSPLRamp, outdev.Fs);

			% plot the stimulus
			subplot(211)
			plot(tonestim);
			drawnow

			% loop through reps for this tone value
			for rep = 1:cal.ToneSPLReps
				% play stimulus, get response;
				resp = arraysingleIO(outdev, indev, zBUS, tonestim, outdev.channel, in_samples);

				% store the data
				speakercal.ToneSPLData{rms_index, tone_index, rep} = resp;

				% filter the data
				respfilt = filter(cal.coeffb, cal.coeffa, resp(tonestart:toneend));
				% compute peak magnitude (in Volts)
				mag = fitsinvec(respfilt, 1, indev.Fs, tonefreq);
				% convert to dB SPL  (take RMS, convert to Pascals, then dB)
				tonedb = dbspl(VtoPa*RMSsin*mag);

				% plot the response
				subplot(212)
				plot(respfilt);
				drawnow

				disp(sprintf('\t \t ...rep %d \t dB SPL = %f', rep, tonedb));

			end	% end of REPS loop
		end	% end of TONES loop
	end	% end of RMS loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% pre-allocate arrays
	mags = zeros(cal.Ntonesplrms, cal.Ntonesplfreqs, cal.ToneSPLReps);
	phis = zeros(cal.Ntonesplrms, cal.Ntonesplfreqs, cal.ToneSPLReps);
	distmags = zeros(cal.Ntonesplrms, cal.Ntonesplfreqs, cal.ToneSPLReps);
	distphis = zeros(cal.Ntonesplrms, cal.Ntonesplfreqs, cal.ToneSPLReps);
	
	% store the tone freq and Vrms vectors
	caldata.tonerms.freqs = cal.ToneSPLFreqs;
	caldata.tonerms.v_rms = cal.ToneSPLRMS;

	for n = 1:cal.Ntonesplrms
		for t = 1:cal.Ntonesplfreqs
			tonefreq = cal.ToneSPLFreqs(t);
			for r = 1:cal.ToneSPLReps
				resp = speakercal.ToneSPLData{n, t, r}(tonestart:toneend);
				respfilt = filter(cal.coeffb, cal.coeffa, resp);

				[mag, phi] = fitsinvec(respfilt, 1, indev.Fs, tonefreq);
				mags(n, t, r) = mag;
				phis(n, t, r) = phi;
				[mag, phi] = fitsinvec(respfilt, 1, indev.Fs, 2*tonefreq);
				distmags(n, t, r) = mag;
				distphis(n, t, r) = phi;
			end	% end of REPS loop
			caldata.tonerms.pa_avg(n, t) = mean(VtoPa*RMSsin*mags(n, t, :));
			caldata.tonerms.pa_std(n, t) = std(VtoPa*RMSsin*mags(n, t, :));
			caldata.tonerms.phi_avg(n, t) = mean(VtoPa*RMSsin*phis(n, t, :));
			caldata.tonerms.phi_std(n, t) = std(VtoPa*RMSsin*phis(n, t, :));
		end	% end of TONES loop
	end	% end of RMS loop

	speakercal.ToneSPLResp.mags = mags;
	speakercal.ToneSPLResp.phis = phis;
	speakercal.ToneSPLResp.distmags = distmags;
	speakercal.ToneSPLResp.distphis = distphis;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure(tonesplfig);

	plot(caldata.tonerms.freqs, dbspl(caldata.tonerms.pa_avg), '.-')

	for n = 1:cal.Ntonesplrms
		legstr{n} = num2str(cal.ToneSPLRMS(n));
	end

	titlestr{1} = sprintf('dB SPL (Tones) vs. Freq (Hz), az: %d el: %d', az, el);
	titlestr{2} = sprintf('Vrms = %.2f - %.2f', ...
							cal.ToneSPLRMS(1), cal.ToneSPLRMS(end));
	title(titlestr)
	xlabel('Tone Frequency (Hz)');
	ylabel('Output dB SPL');
	grid on
	legend(legstr);
	set(gcf, 'Name', sprintf('Az: %d El %d Tone Response Function', az, el));
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_tonespl' '.fig'];
	saveas(tonesplfig, figfilename)

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FFCal_ToneSPLResponse.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
