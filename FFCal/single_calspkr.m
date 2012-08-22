function spkrinfo = single_calspkr(az, el)

spkrinfo = [];

sdata = speaker_data('azel', [az el]);

if sdata.speakernum == -1
	error('speaker not found - exiting!!!!!');
end

spkr = 1;
spkrinfo{spkr} = sdata;
spkrinfo{spkr}.calibrated = 0;
fname = sprintf('%d_%d', az, el);
fname(find(fname=='-')) = 'n';
spkrinfo{spkr}.datafile = fname;
		