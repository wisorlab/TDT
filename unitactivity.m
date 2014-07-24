function output = unitactivity( params )

	
	% find unit activity in data recorded by the TDT system
	% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	% 
	% :param params - struct
	% params is a struct with the following fields:
	% 
	% :field channel - channel to look at
	% :field sort = PCA sortcode(s) referring to a cluster of units to be proccessed in parallel 
	% with the EEG. The sort field can be a scalar or vector of multiple codes.
	% :field snippetWidth - width (in seconds) of the snippet window
	% :field snippetOffset - offset (in seconds) into the file to start processing
	% :field plot - boolean value of whether or not to plot the sort codes
	% 
	% --- If the data field is provided, the fields 'tank' and 'block' are not needed and will be ignored 
	% in favor of the data passed in as 'data' ---
	% 
	% :field data - (optional) If a TDT block data set already exists as a matlab structured array, the
	% function may plot data from that structured array (passed in as 'data'), otherwise the parameters 
	% 'tank' and 'block' must be provided.
	% \
	% --- If the following parameters are provided, this function will generate data from the specified
	% tank and block and (optionally) plot the resulting data. ---
	% 
	% :field tank - path to a TDT tank
	% :field block - block name within 'tank'

	% ~~~~~~~~~~~~~~~~~~~~
	% parse the parameters
	% ~~~~~~~~~~~~~~~~~~~~
	
	defaults = struct( ...
		'channel', 1, ...
		'sort', [1:3], ...
		'snippetWidth',50, ...
		'snippetOffset',50, ...
		'data',[], ...
		'plot','true');
  
	params = extend(defaults,params);
 
	if ~isfield(params,'data') && ~empty(params.data)
		params.data = TDT2mat(tank,block);
	end
	if ~isfield(params,'sort') && ~empty(params.sort)
		params.sort = 1:3;
	end
	
	% make calling these easier
	eNeu = params.data.snips.eNeu;
	EEGx = params.data.streams.EEGx;
	snippetOffset = params.snippetOffset;
	snippetWidth = params.snippetWidth;
	
	for i=1:length(params.sort)
		chansort{i} = eNeu.ts(intersect(find(eNeu.chan==params.channel),find(eNeu.sortcode==params.sort(i))));
	end
	
	if snippetOffset < 1/EEGx.fs
		snippetOffset = 1/EEGx.fs;
	end

	output = struct();
	output.window = round(EEGx.fs * snippetOffset):round(EEGx.fs * snippetOffset + EEGx.fs * snippetWidth);
	output.x = (output.window/EEGx.fs);
	output.y = double(EEGx.data(params.channel,output.window));
	output.fs = EEGx.fs;
	output.data = params.data;

	% if the plot parameter is set to true plot the data
	if ~strcmp(params.plot,'false')

		figure
		set(gca,'Color','black');
		hold on
		plot(output.x,output.y,'w')
		a = [snippetOffset (snippetOffset + snippetWidth) (-max(EEGx.data(output.window))*1.5) max(EEGx.data(output.window))];
		axis(a)
		title(['EEG Signal in channel ' num2str(params.channel) ])

	end

	output.sortEvents = cell(length(params.sort));
		
	% calculate (and optionally plot) sortcodes
	colors = {'r.','m.','g.'};
	for i=1:length(params.sort)
		sortEvents = chansort{:,i}(chansort{:,i} > snippetOffset & chansort{:,i} < snippetOffset + snippetWidth);
		output.sortEvents{i} = sortEvents;

		% if the plot parameter is set to true plot the data
		if ~strcmp(params.plot,'false')
			plot(sortEvents,ones(length(sortEvents)).*a(3)+(abs(a(3)-a(4))*0.05*i),colors{i})
			axis(a)
			title(['Sort Code ' num2str(i)])
		end
	end
end