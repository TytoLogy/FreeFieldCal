%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup data storage variables and
% paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	datadir = date;
	init_spkrinfo = 1;

	tab = sprintf('\t');
	ret = sprintf('\n');
	
	options_list = {	[ret 'FFCal Options:' ret], ...
						[tab '1' tab tab 'Calibrate All Speakers' ], ...
						[tab '2' tab tab 'Calibrate Single Speaker' ], ...
						[tab '3' tab tab 'Abort program' ], ...
						[ret]	};

	for i = 1:length(options_list)
		disp(options_list{i});
	end
	
	q = query_userint('Select option', [1 3 2]);
	
	ABORT = 0;
	
	switch q
		
		case 1
			% check to see if the spkrinfo.mat file exists
			if exist('spkrinfo.mat')
				disp('***Previous spkrinfo.mat exists!***');
				q1 = query_user('Continue at last speaker', 1);
				if q1 == 0
					disp('Warning!')
					disp('Previously collected data will be overwritten!')
					q2 = query_user('Continue with overwrite?', 0);
					if q2
						init_spkrinfo = 1;
					end
				end

				if q1 == 1 
					% if spkrinfo.mat exists and user chooses to continue as 
					% previous point, load the old spkrinfo.mat data
					init_spkrinfo = 0;
					load('spkrinfo.mat');

					% find where we left off
					i = 1;
					while spkrinfo{i}.calibrated
						i = i+1;
					end
					start_speaker = i;
				end
			else
				init_spkrinfo = 1;
			end

			% otherwise (or if user chooses), start from beginning
			if init_spkrinfo
				disp('Creating spkrinfo')
				if exist('spkrinfo.mat')
					delete('spkrinfo.mat');
				end
				spkrinfo = calspkrorder;
				save spkrinfo.mat spkrinfo
				start_speaker = 1;
			end

			% create the data directory
			if exist(datadir, 'dir')
				disp([mfilename ': data directory ' datadir ' exists.']);
			else
				disp([mfilename ': creating data directory ' datadir '...'])
				mkdir(datadir);
			end

	
		case 2
			useraz = query_int('Enter Azimuth of speaker to calibrate', [-100 100]);
			userel = query_int('Enter Elevation of speaker to calibrate', [-80 80]);
			spkrinfo = single_calspkr(useraz, userel);
			start_speaker = 1;

			% create the data directory
			if exist(datadir, 'dir')
				disp([mfilename ': data directory ' datadir ' exists.']);
			else
				disp([mfilename ': creating data directory ' datadir '...'])
				mkdir(datadir);
			end

		case 3
			ABORT = 1;
			return;
			
		otherwise
			return;
	end
	