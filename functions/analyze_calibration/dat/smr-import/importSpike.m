function datout = importSpike(filename, chans)
%% Converts Spike File to continuous data structure
%
% Input:
%   d - data structure returned frmo readSpikeFile.m
%   (chans) - which channels to import
%
% Output:
%   datout - dat object with fields:
%     'data' - column vector of numeric data (row)
%     'chanlabels' - string labels for each channel
%     'chanvals' - numeric values associated with each
%         channel (e.g. frequency at each bin for a spectrogram).
%     'samplerate' - derived, not claimed (samples/second, double-precision)
%     'tstart','tend' - time of first/last sample; in seconds.
%     'units' - descriptive string giving units of data.
%     'nbad_start' and 'nbad_end' - samples at start/end of data that are
%         unreliable due to filtering edge effects)
%
% Dependencies: SON library, readSpikeFile.m
%
% Hannah Payne July 2012
%

%% Load data
if ~exist('chans', 'var')
    chanlist = readSpikeFile(filename,[]);
    chans = [chanlist(:).number];    
end


if iscell(chans)
    chanlist = readSpikeFile(filename,[]);
    temp = arrayfun(@(x) any(strcmp(x.title,chans)),chanlist);
    chans = [chanlist(temp).number];
end

d = readSpikeFile(filename,chans);
nchans = length(d);

% Preallocate object array
datout(nchans) = dat;

%% Load data
iIn=0;
iOut = 0;
itend = 0;

while iIn < nchans
    iIn=iIn+1;
    iOut=iOut+1;
    ichanval = d(iIn).header.channel;
    ichanlabel = d(iIn).header.title;
    
    if strcmp(d(iIn).header.channeltype,'Episodic Waveform')
        d(iIn) = fillGaps(d(iIn));
    end
    
    % Continuous data
    if strcmp(d(iIn).header.channeltype,'Continuous Waveform')
        % Scale current data
        iData = single(d(iIn).data.adc);
        scaleFactor = d(iIn).header.adc.Scale;
        iData = scaleFactor*iData;
        
        % Store it
        isamplerate = 1/prod(d(iIn).header.adc.SampleInterval);
        itstart = 0;
        itend = (length(iData)-1)/isamplerate;
        %     itend = d(i).data.tim(end) - d(i).data.tim(1);
        iunits = d(iIn).header.adc.Units;
        
        % Event date
    else
        iData = d(iIn).data.tim;
        
        if ~isempty(d(iIn).data.mrk) && any( d(iIn).data.mrk(:) ~=0 )
            isamplerate = native2unicode(d(iIn).data.mrk(:,1));
        else
            isamplerate = 'event';
        end        
        iunits = 's';
        itstart = 0;
        itend = max(itend,max(iData));
    end
    
    % Make the object instance
    if size(iData,2)==2
        datout(iOut) = dat(iData(:,1),ichanlabel,ichanval,isamplerate,itstart,itend,iunits);
        datout(iOut+1) = dat(iData(:,2),ichanlabel,ichanval,isamplerate,itstart,itend,iunits);
        if iscell(datout(iOut).chanlabel)
            currlabel = datout(iOut).chanlabel{:};
        else
            currlabel = datout(iOut).chanlabel;
        end        
        datout(iOut).chanlabel = [currlabel '+'];        % Rising edge
        datout(iOut+1).chanlabel = [currlabel '-'];    % Falling edge
        iOut = iOut+1;
    else
        datout(iOut) = dat(iData,ichanlabel,ichanval,isamplerate,itstart,itend,iunits);
    end
end

%% Make sure all channels with same samplerate are the same size
contchans = arrayfun(@(x) isnumeric(x.samplerate),datout);
datoutcont = datout(contchans);

samplerates = arrayfun(@(x) x.samplerate, datoutcont);
ns = arrayfun(@(x) length(x.data),datoutcont);
mask = samplerates == mode(samplerates);
targetN = min(ns(mask));
for i = 1:length(datoutcont)
    if mask(i)
    datoutcont(i).data = datoutcont(i).data(1:targetN);
    end
end

datout(contchans) = datoutcont;

end




function dout = fillGaps(din)
%% dout = fillGaps(din) Fill in gaps in older Spike2 recordings
% Hannah Payne 12/16/13 

dout = din;
nsegs = size(din.data.tim,1);
maxt = din.data.tim(end); % seconds
samplerate = 1/prod(din.header.adc.SampleInterval);
npoints = round(maxt*samplerate+1);
data = zeros(1,npoints);
tt = (1:npoints)/samplerate;

for j = 1:nsegs
    currlength = din.header.adc.Npoints(j);
    currdata = din.data.adc(1:currlength,j);
    currtime = din.data.tim(j,1);
    startind = find(tt>=currtime,1);
    data(startind:startind+currlength-1) = currdata;
end

dout.data.adc = data;
dout.data.tim = [0 maxt];
dout.header.channeltype = 'Continuous Waveform';
end




