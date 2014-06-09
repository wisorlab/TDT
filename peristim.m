bin_size = 0.05; % (in seconds)

data = TDT2mat('C:\Users\wisorlab\Desktop\6_9_14_TEST','5min_Detect_Units');


ttl_ts = data.epocs.Stim.onset;
eNeu = data.snips.eNeu.chan;
eNeu_ts = data.snips.eNeu.ts;

% add up the number of events in each channel in "chans" and "randChans"
chans = zeros(16,1);
randChans = zeros(size(chans));

% create a random matrix the same size
randoms = floor( rand( size(eNeu) ) * 16 + 1 );


for i=1:length(ttl_ts)

    % select a duration of "bin_size" after each stimulus onset
    range = eNeu( eNeu_ts > ttl_ts(i) & eNeu_ts < (ttl_ts(i)+ bin_size));
    randRange = randoms( eNeu_ts > ttl_ts(i) & eNeu_ts < (ttl_ts(i)+ bin_size));
    
    %increment chans by the number of events in each channel
    for j=1:length(chans)
        randChans(j) = randChans(j) + length(randRange( randRange == j ));
        chans(j) = chans(j) + length(range( range == j ));
    end
    
end


figure
hold on
bar(1:length(chans),chans,1,'r');
bar(1:length(randChans),randChans,0.5,'b')

xlabel('Channel')
ylabel('Number of Events')
legend('eNeu','Random')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3D Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% figure
% bar3([ chans(1:8), chans(9:end)])
