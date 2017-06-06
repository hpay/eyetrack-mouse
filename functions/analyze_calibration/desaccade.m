function eyevelpost = desaccade(eyevel, samplerate, varargin)
% EYEVELPOST = desaccade(eyevel, samplerate, threshold, ploton,'window',window)

p = inputParser;
addOptional(p,'threshold',50);
addOptional(p,'ploton',0);
addParamValue(p,'window',[.1 .3]);
addParamValue(p,'freq',[]);

parse(p,varargin{:})
ploton = p.Results.ploton;
threshold = p.Results.threshold;
window = p.Results.window;
freq = p.Results.freq;
% time = p.Results.time;
eyevel = double(eyevel(:));
time = (1:length(eyevel))'/samplerate;

if isempty(threshold)
    threshold = 50;
end

window = abs(window);
presaccade = ceil(window(1)*samplerate); % Fixed bug - divide by sr 7/21/13
postsaccade = ceil(window(2)* samplerate); 

eyevel2 = eyevel(:);
eyevel2(abs(eyevel2)>threshold*2)=NaN;
if ~isempty(freq)
[amp, phase] = fitsine(time, eyevel2, freq,0);
eyevel3 = eyevel-amp*sin(time*2*pi*freq+deg2rad(phase));
else
    eyevel3 = eyevel;
end

%% Find saccades in eye movement and blank out an interval on either side
%find regions of data outside threshold; 
rejecttemp0 = abs(eyevel3) > threshold | isnan(eyevel3);

%remove points around omit centers as defined by pre- & postsaccade time
sacmask = ones(1,presaccade+postsaccade);

%filter function replaces zeros with ones (equal to remove time) around an omit center
rejecttemp1 = conv(double(rejecttemp0),sacmask,'full');
rejecttemp2 = rejecttemp1(presaccade:presaccade+length(eyevel3)-1);

% eyevel with desaccade segments removed
eyevelpost = eyevel;
eyevelpost(logical(rejecttemp2))= NaN;

%% DEBUG
if ploton
%     figure(ploton);clf;

    plot(time,eyevel,'k','LineWidth',1); hold on
    plot(time, eyevelpost,'r','LineWidth',.5);
%     try
%         plot(find(rejecttemp0)/samplerate,0,'ob','LineWidth',3)
%     catch
%     end
    ylim([-50 50])    
end
