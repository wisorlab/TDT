
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% detect unit activity and slow wave patterns
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

addpath './tdt'
tank = 'D:\OptoTagging\7_23_14_R2180';
block = 'Stimulate_10min';

data = TDT2mat(tank,block,'VERBOSE',false);

chan = 1;
sortCode = 2;
snippetWidth = 100;
snippetOffset = 100;

% get the unit activity
units = unitactivity( struct( ...
	'channel', chan, ...
	'sort', sortCode, ...
	'snippetWidth',snippetWidth, ...
	'snippetOffset',snippetOffset, ...
	'data',data, ...
	'plot', 'false' ));


% set up parameters that define a slow wave:
% 
% - Chebyshev Type II Parameters
% 	- Passband Edges
% 		- lowPassEdge (?)
% 		- highPassEdge (?)
% 	- StopBand Edges
% 		- lowStopEdge (?)
% 		- highStopEdge (?)
% 	- passBandRipple (?)
% 	- stopBandAttenuation (?)
% 
chb = struct( ...
	'lowPassEdge', 0.5, ...
	'highPassEdge', 4, ...
	'lowStopEdge', 0.01, ...
	'highStopEdge', 10, ...
	'passBandRipple', 3, ...
	'stopBandAttenuation', 20 );

% - Wave Detection Parameters
% 
% 	- amplitude (microV)
% 	- maxAmplitude (microV)
% 	- p2ptime (seconds)
% 	- maxp2ptime (seconds)
% 	- epochLength (seconds)
% 
params = struct( ...
	'amplitude', 0.100, ...
	'maxAmplitude', 0.600, ...
	'p2ptime', 0.25, ...
	'maxp2ptime', 2, ...
	'epochLength', 4 );

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% slow wave detection algorithms
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Chebyshev filtering
Wp=[ chb.lowPassEdge  chb.highPassEdge ]/(units.fs/2); Ws=[ chb.lowStopEdge chb.highStopEdge ]/(units.fs/2);  
[n, Wn]=cheb2ord(Wp,Ws,chb.passBandRipple, chb.stopBandAttenuation);
[bb,aa]=cheby2( n, chb.stopBandAttenuation, Wn );

try
	filtered_matrix = filtfilt(bb,aa,units.y);
catch err
	warning(['The current parameters create a non-invertible matrix that would create improper scaling' ...
		'Vary the Chebyshev parameters or compression.']);
end


% the 6th column represents the voltage difference between the maximum and the first point in the wave.
slowWaves = detectSlowWaves( filtered_matrix, units.fs, params.epochLength );
slowWaves_ = slowWaves( slowWaves(:,5)>=params.p2ptime,:);
slowWaves_ = slowWaves( slowWaves(:,5)<=params.maxp2ptime,:);
slowWaves_ = slowWaves( slowWaves(:,6)>=params.amplitude,:);
slowWaves_ = slowWaves( slowWaves(:,6)<=params.maxAmplitude,:);
tableWaves = slowWaves;
% tableWaves(:,2:4) = tableWaves(:,2:4)/units.fs;



