
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% detect unit activity and slow wave patterns
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% p1=str2double(get(handles.lowPassEdge,'String')); 
% p2=str2double(get(handles.highPassEdge,'String')); % pass band
% s1=str2double(get(handles.lowStopEdge,'String'));
% s2=str2double(get(handles.highStopEdge,'String')); % stop band
% Rp=str2double(get(handles.passBandRipple,'String')); 
% Rs=str2double(get(handles.stopBandAttenuation,'String'));


% Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2);  
% [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
% [bb,aa]=cheby2(n,Rs,Wn);
% filtered_matrix = handles.matrix;
% s = warning('error','MATLAB:nearlySingularMatrix');
% try
%     filtered_matrix(:,2:end) = filtfilt(bb,aa,handles.matrix(:,2:end));
% catch err
%     msgbox('The current parameters create a non-invertible matrix that would create improper scaling. Vary the Chebyshev parameters or compression.');
% end
% warning(s);

% tdt.TDT2Mat('data','6_17_14_10min_Stim');

tank = 'data/6_17_14_10min_Stim';
sev = '6_17_14_6_17_14_10min_Stim_EEGx_Ch'
chan = 1;

sevData = tdt.SEV2mat(tank,'CHANNEL',chan);

eeg = sevData.EEGx;

detectSlowWaves(eeg.data,eeg.fs,10);