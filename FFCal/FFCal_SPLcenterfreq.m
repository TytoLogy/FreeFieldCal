% FFCal_SPLcenterfreq.m
% 
% Tests how SPL varies for narrowband noise as a function of center 
% frequency of calibrated stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created:
%	8 June, 2009
%
% Revisions:
%------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the FFCal_SPLcenterfreq Curve Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Running SPL - Fctr determination');
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
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for SPL NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	arraySetTiming(outdev, indev, cal.SPLSweepPeriod, cal.SPLDelay,...
									cal.SPLDuration, cal.SPLAcquisition);
	% Set the sweep count (may not be necessary)
	RPsettag(outdev, 'SwCount', 1);
	RPsettag(indev, 'SwCount', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play narrowband noise at different Fctr (if specified) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	nBW = length(cal.BWnarrow.fctr);
	if nBW >= 1
		disp('dB vs Fctr for Narrowband stimuli');

		BWctrspl = zeros(nBW, cal.SPL_BWreps);

		for n = 1:nBW
			fctr = cal.BWnarrow.fctr(n);
			bw = cal.BWnarrow.bw(1);
			
			disp(sprintf('...BW = %d \t Fctr = %d', bw, fctr));
			
			% synthesize the broadband stimulus
			fmin = fctr - bw / 2;
			fmax = fctr + bw / 2;

			% get the stimulus intensity level scaling factor from the 
			% calibration data
			% NOTE that we need to do this at each bw, because it is
			% bw dependent!!!!!
			scale = get_bw_scale_nowarning(cal.BWtestlevel, caldata, [fmin fmax]);
			disp(['bw scale factor = ' num2str(scale)]);

			[stim, mag, phase]  = synmononoise_fft(cal.SPLDuration, outdev.Fs, fmin, fmax, scale, caldata);		
			% ramp the stimulus on and off for smoothness' sake...
			stim = sin2array(stim, cal.TestRamp, outdev.Fs);

			for r = 1:cal.SPL_BWreps
				% play the test narrowband stimulus
				[resp, index] = arraysingleIO(outdev, indev, zBUS, stim, outdev.channel, in_samples);

				%filter the test stimulus (eliminate the nasty low frequency crap)
				testresp = filter(cal.coeffb, cal.coeffa, resp);
				testresp = testresp(cal.splstart:cal.splend);
				testdbspl = dbspl(rms(VtoPa*testresp));

				disp(sprintf('\t\t rep %d \t Desired dB SPL: %.2f \t Measured dB SPL: %.2f',...
									r, cal.BWtestlevel, testdbspl));
				
				BWctrspl(n, r) = testdbspl;

				fftplot(testresp, indev.Fs, stimrespfig1);
				subplot(311)
				title('Frequency Response, NB CTR');
			end
		end

		% compute mean
		caldata.BWctrfreqs = cal.BWnarrow.fctr;
		caldata.BWctrdB = mean(BWctrspl, 2);
		caldata.BWctrdB_std = std(BWctrspl, 0, 2);
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% plot
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		figure(fctrsplfig)
		errorbar(caldata.BWctrfreqs, caldata.BWctrdB, caldata.BWctrdB_std, 'b.-')
		xlabel('Fctr (Hz)')
		ylabel('dB SPL');
		titlestr{1} = sprintf('dB vs. Fctr, Narrowband Noise, az: %d el: %d', az, el);
		titlestr{2} = sprintf('BW = %d, Target dB = %d', ...
							cal.BWnarrow.bw(1), cal.BWtestlevel);
		title(titlestr)
		grid on
		set(gcf, 'Name', 'dB vs. Fctr, NB');
		figfilename = [datadir '\' spkrinfo{spkr}.datafile '_xferfunction' '.fig'];
		saveas(fctrsplfig, figfilename)

	end % END of variable ctr freq IF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FFCal_SPLcenterfreq.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
