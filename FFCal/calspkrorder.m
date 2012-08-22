function spkrinfo = calspkrorder;

spkr = 0;
start_speaker = 1;
for az=-100:10:100	
	% look for the els at given azimuth
	sdata = speaker_data('azimuth', az);
	Nels = length(sdata.elevations);
	if Nels
		% sort the elevations from top to bottom
		els = sort(sdata.elevations, 'descend');
		for el = 1:Nels
			% increment speaker counter;
			spkr = spkr + 1;

			spkrinfo{spkr} = speaker_data('azel', [az els(el)]);
			spkrinfo{spkr}.calibrated = 0;
			fname = sprintf('%d_%d', az, els(el));
			fname(find(fname=='-')) = 'n';
			spkrinfo{spkr}.datafile = fname;
		end
	end
end
		