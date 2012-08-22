%------------------------------------------------------------------------
% FFCal_ToneResponse.m
%------------------------------------------------------------------------
%
% Using calibration from B&K mic
% Tests tone generation
%
%------------------------------------------------------------------------
% Part of FFCal program
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created:
%	16 April, 2009
%
% Revisions:
%	29 May, 2009 (SJS):	finished conversion from ClickResponse code,
%									testing to begin
%	8 June, 2009 (SJS):	recreated
%	28 October, 2009 (SJS): more documentation
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play tones to check for proper calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Running Tone Response...');
	pause(0.5);
	
	% compute # of samples to collect
	in_samples = ms2samples(cal.ToneAcquisition, indev.Fs);
	% start bin for tone signal
	tonestart = ms2bin(cal.ToneDelay + cal.ToneRamp, indev.Fs);
	% end bin for tone signal
	toneend = tonestart + ms2bin(cal.ToneDuration, indev.Fs);
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for TONES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	status = arraySetTiming(outdev, indev, cal.ToneSweepPeriod, cal.ToneDelay,...
									cal.ToneDuration, cal.ToneAcquisition);
	% Set the sweep count (may not be necessary)
	RPsettag(outdev, 'SwCount', 1);
	RPsettag(indev, 'SwCount', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop, play tones, record resps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create figure
	figure(stimrespfig1);
	
	% loop through frequencies
	for tone_index = 1:cal.Ntonefreqs
		% get the tone frequency to test
		tonefreq = cal.ToneFreqs(tone_index);
		
		disp(sprintf('...tone %.2f response, %d of %d', tonefreq, tone_index, cal.Ntonefreqs));
		pause(0.1);

		% get the scale factor for this frequency
		tonescale = get_bw_scale_nowarning(cal.ToneLevel, caldata, tonefreq);
		disp(sprintf('\t \t tonelevel = %d tonescale = %.4f', cal.ToneLevel, tonescale));
		% store in speakercal structure element
		speakercal.ToneResp.scale(tone_index) = tonescale;
		
		% synthesize, window and plot the tone
		tonestim = synmonosine(cal.ToneDuration, outdev.Fs, tonefreq, tonescale, caldata);
		tonestim = sin2array(tonestim, cal.ToneRamp, outdev.Fs);
		subplot(211), plot(tonestim);
		
		% loop through desired # reps
		for n = 1:cal.ToneReps
			% play stimulus, get response;
			[toneresp, index] = arraysingleIO(outdev, indev, zBUS,...
												tonestim, outdev.channel, in_samples);
			% filter response
			tonerespfilt = filter(cal.coeffb, cal.coeffa, toneresp(tonestart:toneend));
			% plot response
			subplot(212)
			plot(filter(cal.coeffb, cal.coeffa, toneresp));
			drawnow
			
			% compute magnitudes and phases for tone, store in structure
			[tmp1, tmp2] = fitsinvec(tonerespfilt, 1, indev.Fs, tonefreq);
			speakercal.ToneResp.mag(tone_index, n) = tmp1;
			speakercal.ToneResp.phase(tone_index, n) = tmp2;
			
			% compute magnitudes and phases for 2nd harmonic distortion, store in structure			
			[tmp1, tmp2] = fitsinvec(tonerespfilt, 1, indev.Fs, 2*tonefreq);
			speakercal.ToneResp.distmag(tone_index, n) = tmp1;
			speakercal.ToneResp.distphase(tone_index, n) = tmp2;
			
			% compute db level of tone
			tonedb = dbspl(VtoPa.*RMSsin.*speakercal.ToneResp.mag(tone_index, n));
			disp(sprintf('\t \t dB SPL = %f', tonedb));			
		end	% end of N (reps) loop
	end	% end of tone_index (freq) loop
	
	% new figure
	figure
	
	% compute means and std devs of measures, store in the caldata struct
	if cal.ToneReps > 1
		caldata.tone.mag_avg = mean(dbspl(VtoPa.*RMSsin.*speakercal.ToneResp.mag'));
		caldata.tone.mag_std = std(dbspl(VtoPa.*RMSsin.*speakercal.ToneResp.mag'));
		tonerror_avg =  mean(dbspl(VtoPa.*RMSsin.*speakercal.ToneResp.mag') - cal.ToneLevel);
		tonerror_std =  std(dbspl(VtoPa.*RMSsin.*speakercal.ToneResp.mag') - cal.ToneLevel);
	else
		caldata.tone.mag_avg = dbspl(VtoPa.*RMSsin.*speakercal.ToneResp.mag');
		caldata.tone.mag_std = zeros(size(caldata.tone.mag_avg));
		tonerror_avg =  dbspl(VtoPa.*RMSsin.*speakercal.ToneResp.mag') - cal.ToneLevel;
		tonerror_std =  zeros(size(tonerror_avg));
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure(tonerespfig)
	subplot(211)
	errorbar(cal.ToneFreqs, caldata.tone.mag_avg, caldata.tone.mag_std, '.-');
	ylabel('dB SPL')
	ylim([0 60])
	titlestr{1} = sprintf('dB vs. Tone Freq, az: %d el: %d', az, el);
	titlestr{2} = sprintf('Target dB = %d', cal.ToneLevel);
	title(titlestr)
	grid on

	subplot(212)
	errorbar(cal.ToneFreqs, tonerror_avg, tonerror_std, '.-');
	grid on
	xlabel('Tone Freq (Hz)')
	ylabel('Error dB SPL')
	ylim([-10 10])
	
	set(gcf, 'Name', 'Calibrated Tone Response');
	figfilename = [datadir '\' spkrinfo{spkr}.datafile '_caltoneresp' '.fig'];
	saveas(tonerespfig, figfilename)
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrap it up... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('Tone Response Finished for speaker [%d, %d]', az, el))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FFCal_ToneResponse.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

