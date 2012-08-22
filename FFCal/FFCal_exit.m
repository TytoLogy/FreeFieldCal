% FFCal_exit.m
%------------------------------------------------------------------------
% script that closes TDT devices and GS3 Multiplexer
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sharad.shanbhag@einstein.yu.edu
%------------------------------------------------------------------------
% Created: ????
%
% Revisions:
%	18 March, 2010 (SJS):
% 		-	made changes to comply with modifications made to some of the RP
% 			toolbox functions
%		-	updated documentation/comments
% 		-	now explicitly closes zBUS device
%------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up the RP circuits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('...closing TDT devices...');

if exist('outdev', 'var')
	disp('...outdev');
	status = RPclose(outdev);
end
if exist('indev', 'var')
	disp('...indev');
	status = RPclose(indev);
end
if exist('servo', 'var')
	disp('...servo');
	status = RPclose(servo);
end
if exist('zBUS', 'var')
	disp('...zBUS');
	status = zBUSclose(zBUS);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up the multiplexer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('...closing GS3 Multiplexer...');
s1 = instrfind;
if exist('gs3obj', 'var')
	try
		gs3_close(gs3obj)
		delete(gs3obj)
	catch
		disp('forcing gs3 close')
		t1 = gs3_readstatus(gs3obj);
		save(['gs3status_' date '.mat'], 't1', 's1');
		gs3_forceclose
	end
end

