if plot_init
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Set up the plots
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	stim_t = 1000 * [0:length(cal.A)-1] ./ outdev.Fs;
	resp_t = 1000 * [0:length(ar)-1] ./ indev.Fs;
	figure
	
	subplot(4, 1, 1)
	stim_a = plot(stim_t, cal.A, 'XDataSource', 'stim_t', 'YDataSource', 'A');
	title('Golay Calibration');
	ylabel('Stim A (V)');
	subplot(4, 1, 2)
	resp_a = plot(resp_t, ar, 'XDataSource', 'resp_t', 'YDataSource', 'ar');
	ylabel('Resp A (V)')
	subplot(4, 1, 3)
	stim_b = plot(stim_t, cal.B, 'XDataSource', 'stim_t', 'YDataSource', 'B');
	ylabel('Stim B (V)');
	subplot(4, 1, 4)
	resp_b = plot(resp_t, br, 'XDataSource', 'resp_t', 'YDataSource', 'br');
	xlabel('Time (msec)');
	ylabel('Resp B (V)');
	
	plot_init = 0;
	
else
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% refresh plots
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	stim_t = 1000 * [0:length(cal.A)-1] ./ outdev.Fs;
	resp_t = 1000 * [0:length(ar)-1] ./ indev.Fs;
	refreshdata(stim_a);
	refreshdata(resp_a);
	refreshdata(stim_b);
	refreshdata(resp_b);
	drawnow;
end

