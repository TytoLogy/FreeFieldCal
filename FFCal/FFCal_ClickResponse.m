% FFCal_ClickResponse.m
%
% Plays clicks, records responses. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created:
%	30 November, 2007
%
% Revisions:
% 	8 June, 2009: recreated
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play clicks to check for echoes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Running Click Response...');
	pause(0.5);
	
	% create the click calibration stimulus
	cal.click = cal.ClickAmplitude*syn_click(cal.ClickDuration, cal.ClickStimDelay, outdev.Fs)';

	% compute # of samples to collect
	in_samples = ms2samples(cal.ClickAcquisition, indev.Fs);
	% start bin for click signal
	clickstart = ms2bin(cal.ClickDelay, indev.Fs);
	% end bin for click signal
	clickend = clickstart + ms2bin(cal.ClickDuration, indev.Fs);
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TDT hardware timing variables for CLICKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	arraySetTiming(outdev, indev, cal.ClickSweepPeriod, cal.ClickDelay,...
									cal.ClickDuration, cal.ClickAcquisition);
	% Set the sweep count (may not be necessary)
	RPsettag(outdev, 'SwCount', 1);
	RPsettag(indev, 'SwCount', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play clicks to check for echoes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for click_index = 1:cal.Nclicks
		disp(sprintf('...click response, %d of %d', click_index, cal.Nclicks));
		pause(1);
		% play stimulus;
		[clkresp, index] = arraysingleIO(outdev, indev, zBUS, cal.click, outdev.channel, in_samples);

		speakercal.Click{click_index} = clkresp;

		clkstim_t = 1000 * (0:length(cal.click)-1) ./ outdev.Fs;
		clkresp_t = 1000 * (0:length(clkresp)-1) ./ indev.Fs;

		if click_index == 1
			% Set up the plots
			figure(stimrespfig1);
			subplot(2, 1, 1)
			clkstimplot = plot(clkstim_t, cal.click, 'XDataSource', 'clkstim_t', 'YDataSource', 'cal.click');
			title('Click Calibration');	ylabel('Click Stim (V)');
			xlim([-10 max(clkstim_t)]);
			subplot(2, 1, 2)
			clkrespplot = plot(clkresp_t, clkresp, 'XDataSource', 'clkresp_t', 'YDataSource', 'clkresp');
			xlabel('Time (msec)');		ylabel('Resp (V)');
			clkplot_init = 0;

		else
			figure(stimrespfig1);
			% refresh plots
			refreshdata(clkstimplot);
			refreshdata(clkrespplot);
			drawnow;
		end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrap it up... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	disp(sprintf('Click Response Finished for speaker [%d, %d]', az, el))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FFCal_ClickResponse.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

