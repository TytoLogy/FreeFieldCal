% FFCal_tdtinit.m
%------------------------------------------------------------------------
% script that sets up TDT parameters for the FFCal free-field speaker
% calibration program
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
%		-	using RPrun to start circuits instead of invoke() routine
% 		-	more complete opening method for zBUS
%------------------------------------------------------------------------


disp('...starting TDT hardware...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize zBus control - this is so we can 
% synchronize everything by sending triggers to all 
% units via zBusTrigA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize zBus control
disp('...starting zBUS...')
zBUS.C = [];
zBUS.status = 0;
tmpdev = zBUSinit('GB');
zBUS.C = tmpdev.C;
zBUS.handle = tmpdev.handle;
zBUS.status = tmpdev.status;
%Set zBUS trigger A to LO 
zBUStrigA(zBUS, 0, 0, 6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the TDT devices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize RX8 device 1
tmpdev = RX8init('GB', outdev.Dnum);
outdev.C = tmpdev.C;
outdev.handle = tmpdev.handle;	
outdev.status = tmpdev.status;

% Initialize RX8 device 2
tmpdev = RX8init('GB', indev.Dnum);
indev.C = tmpdev.C;
indev.handle = tmpdev.handle;	
indev.status = tmpdev.status;

% Initialize RX5
tmpdev = RX5init('GB');
servo.C = tmpdev.C;
servo.handle = tmpdev.handle;	
servo.status = tmpdev.status;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads circuits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdev.rploadstatus = RPload(outdev);
indev.rploadstatus = RPload(indev);
servo.rploadstatus = RPload(servo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starts Circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdev.status = RPrun(outdev);
indev.status = RPrun(indev);
servo.status = RPrun(servo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Status1 = RPcheckstatus(outdev);
Status2 = RPcheckstatus(indev);
Status3 = RPcheckstatus(servo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Query the sample rate from the circuit and set up the time vector and 
% stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdev.Fs = RPsamplefreq(outdev);
indev.Fs = RPsamplefreq(indev);
servo.Fs = RPsamplefreq(servo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up some of the buffer/stimulus parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size of the Serial buffer
npts=150000;  

%specify filter function impulse response
dt = 1/indev.Fs;
GolayDuration = 1000 * (length(Astim) / outdev.Fs);
%number of points to write to buffer
obufpts = ms2samples(GolayDuration, outdev.Fs); 
t = 0:dt:(0.001*GolayDuration - dt);

mclock=RPgettag(indev, 'mClock');

% Set the total sweep period time
RPsettag(outdev, 'SwPeriod', ms2samples(cal.SweepPeriod, outdev.Fs));
RPsettag(indev, 'SwPeriod', ms2samples(cal.SweepPeriod, indev.Fs));
% Set the sweep count (may not be necessary)
RPsettag(outdev, 'SwCount', 1);
RPsettag(indev, 'SwCount', 1);
% Set the Stimulus Delay
RPsettag(outdev, 'StimDelay', ms2samples(cal.StimDelay, outdev.Fs));
% Set the Stimulus Duration
RPsettag(outdev, 'StimDur', ms2samples(cal.StimDuration, outdev.Fs));
% Set the length of time to acquire data
RPsettag(indev, 'AcqDur', ms2samples(cal.AcqDuration, indev.Fs));
% set the ttl pulse duration
% RPsettag(outdev, 'PulseDur', ms2samples(cal.TTLPulseDur, outdev.Fs));

% Set the input and output channels
RPsettag(outdev, 'OutputChannel', outdev.channel);
RPsettag(indev, 'InputChannel', indev.channel);

oOutputChannel = RPgettag(outdev, 'oOutputChannel')
oInputChannel = RPgettag(indev, 'oInputChannel')

%Setup filtering if desired
if cal.InputFilter
	RPsettag(indev, 'HPFreq', cal.InputFc);
end


