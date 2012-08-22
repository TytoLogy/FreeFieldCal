%------------------------------------------------------------------------
% SpeakerStats.m
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbha@aecom.yu.edu
%------------------------------------------------------------------------
% Created:	9 June, 2009
% 
% Revisions:
%------------------------------------------------------------------------

warning off MATLAB:log:logOfZero

clear all
close all

% indatadir = 'W:\Sharad\Calibration\FreeField\CalibrationData\ProcessedData\13-Jan-2008\NetCal';
% indatafile = 'netcal-13-Jan-2008.mat';
% load([indatadir '\' indatafile]);

indatadir = 'W:\Sharad\Calibration\FreeField\CalibrationData\ProcessedData\09-Jun-2009\NetCal-09-Jun-2009';
indatafile = 'netcal-09-Jun-2009.mat';
load([indatadir '\' indatafile]);

NSPEAKERS = 144;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, process the Golay Calibration Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Computing Average Transfer function');

	for s=1:NSPEAKERS		
		sdata = netcal{s};
		
		if s==1
			% allocate arrays
			freqs = sdata.freq;
			nfreqs = length(freqs);
			vrms = sdata.v_rms;
			nvrms = length(vrms);
			netmag = zeros(NSPEAKERS, nfreqs);
			netphase = netmag;
			netmaginv = netmag;
			netparms = zeros(NSPEAKERS, nvrms);
			netdbspl = netparms;
		end
		netmag(s, :) = sdata.mag;
		netphase(s, :) = sdata.phase;
		netmaginb(s, :) = sdata.maginv;
		netparms(s, :) = sdata.pa_rms;
		netdbspl(s, :) = sdata.dbspl;
		
		figure(1)
		subplot(221)
		semilogy(freqs, netmag(s, :));
		title(sprintf('Az: %d El: %d',...
			netcal{s}.speakerinfo.azimuths,...
			netcal{s}.speakerinfo.elevations))

		subplot(223)
		plot(freqs, netphase(s, :))
				
		subplot(222)
		plot(vrms, netparms(s, :));
		subplot(224)
		plot(vrms, netdbspl(s, :));
		
		pause
	end

% plot mean response

	magavg = mean(netmag);
	magstd = std(netmag);
	
	
	
	mag_grandavg = mean(mean(db(netmag)))
	mag_grandstd = mean(std(db(netmag)))
	
	
	phaseavg = 1000*mean(netphase)./freqs;
	phasestd = 1000*std(netphase)./freqs;
	
	phase_grandavg = mean(1000*mean(netphase)./freqs)
	phase_grandstd = mean(1000*std(netphase)./freqs)
	
	figure(2)
		subplot(211)
		plot(freqs, db(magavg), 'k-')
		hold on
			plot(freqs, db(magavg + magstd), 'r:')
			plot(freqs, db(magavg - magstd), 'b:')
		hold off
		ylabel('Magnitude (dB)')
		title('Speaker Frequency Response Mean (+/- 1 stdev)')
		grid
		subplot(212)
		plot(freqs, db(phaseavg), 'k-')
		hold on
			plot(freqs, db(phaseavg + phasestd), 'r:')
			plot(freqs, db(phaseavg - phasestd), 'b:')
		hold off
		grid
		xlabel('Frequency (Hz)')
		ylabel('Phase (usec)')
	
	dbavg = mean(netdbspl);
	dbstd = std(netdbspl);
	
	figure(3)
	plot(vrms, dbavg, 'k-')
	hold on
		plot(vrms, dbavg + dbstd, 'r:')
		plot(vrms, dbavg - dbstd, 'b:')
	hold off
	grid

	figure(4)
	plot(freqs, db(netmag), ':')
	hold on
		plot(freqs, db(magavg), 'k-', 'LineWidth', 2)
		plot(freqs, db(magavg + magstd), 'r-', 'LineWidth', 2)
		plot(freqs, db(magavg - magstd), 'b-', 'LineWidth', 2)
	hold off
	
	figure(5)
	plot(vrms, netdbspl, ':', 'LineWidth', 0.1)
	hold on
		plot(vrms, dbavg, 'k-', 'LineWidth', 2)
		plot(vrms, dbavg + dbstd, 'r-', 'LineWidth', 2)
		plot(vrms, dbavg - dbstd, 'b-', 'LineWidth', 2)
	hold off
	grid


warning on MATLAB:log:logOfZero


	