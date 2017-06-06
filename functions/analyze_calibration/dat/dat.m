classdef dat
    % DAT continuous-time AND discrete event data storage object
    %
    %   PROPERTIES
    %     'data' - matrix of numeric data, 1 row per sample, 1 column per
    %         channel. Can be measurement values (continuous) or event times
    %         (discrete)
    %     'chanlabel' - string labels for each channel
    %     'chanval' - numeric values associated with each
    %         channel (e.g. frequency at each bin for a spectrogram).
    %     'samplerate' - samples/second (continuous), 'event' for event channel
    %     'tstart','tend' - time of first/last sample; in seconds.
    %     'units' - descriptive string giving units of data.
    %     'nbadstart' and 'nbadend' - samples at start/end of data that are
    %         unreliable due to filtering edge effects)
    %
    %   METHODS
    %         obj = dat(data, chanlabel, chanval, samplerate, tstart, tend, units)
    %         chanind = datchanind(obj,chans)
    %         datout = datchan(obj, chans)
    %         data = datchandata(obj,chan)
    %         datout = datseg(obj,seg)
    %         data = datsegdata(obj, tstart, duration)
    %         datout = datsegmeans(obj, tstart, duration, shiftsegs,nonan)
    %         datout = dateventtocont(obj, npoints, tstart, tend)
    %         datrange = datrange(obj,chan)
    %         datout = datsmooth(obj,varargin)
    %         datout = datbin(obj, binwidth)
    %         datout = downsample(obj, varargin)         
    %         datout = upsample(obj, n)
    %         datout = resettime(obj);
    %         datout = dattimeind(obj,t) % returns data point at time t
    %         datout = deriv(obj);
    %         t = dattime(obj, seg)
    %         h = plot(obj, varargin)
    %               datfft
    %
    %   DEPENDENCIES
    %       tight_subplot.m
    %
    % Author: Hannah Payne 2012
    %
    properties
        data;
        chanlabel;
        chanval;
        samplerate;
        tstart;
        tend;
        units;
        nbadstart=0;
        nbadend=0;
    end
    methods
        % Class constructor
        function obj = dat(varargin)
            %   dat(data, chanlabel, chanval, samplerate, tstart, tend,
            %   units)
            %     'data' - matrix of numeric data, 1 row per sample, 1 column per
            %         channel.
            %     'chanlabel' - string labels for each channel
            %     'chanval' - numeric values associated with each
            %         channel (e.g. frequency at each bin for a spectrogram).
            %     'samplerate' - derived, not claimed (samples/second, double-precision)
            %     'tstart','tend' - time of first/last sample; in seconds.
            %     'units' - descriptive string giving units of data.
            %
            %     'nbadstart' and 'nbadend' - (not set by user) samples at start/end of data that are
            %         unreliable due to filtering edge effects)
            %
            %    May be called with any or no input arguments
            %
            
            p = inputParser;
            addOptional(p,'data',[],@(x)isnumeric(x) || islogical(x));
            addOptional(p,'chanlabel',{},@(x) ischar(x) || iscell(x) ||isempty(x));
            addOptional(p,'chanval',[],@isnumeric);
            addOptional(p,'samplerate',1,@(x) isnumeric(x) || strcmp(x,'event') || ischar(x));
            addOptional(p,'tstart',0,@isnumeric);
            addOptional(p,'tend',[],@isnumeric);
            addOptional(p,'units',{},@(x) ischar(x) || iscell(x) ||isempty(x));
            
            parse(p,varargin{:});
            
            data_in = p.Results.data;
            chanlabels = p.Results.chanlabel;
            chanvals = p.Results.chanval;
            samplerates = p.Results.samplerate;
            tstarts = p.Results.tstart;
            tends = p.Results.tend;
            unitss = p.Results.units;
            
            if nargin>0
                if size(data_in,1)<size(data_in,2)
                    data_in = data_in';
                end
                if isempty(data_in)
                    nchans = 1;
                else
                    nchans = size(data_in,2);
                end
                
                % Preallocate object array
                obj(nchans) = dat;
                
                if ischar(chanlabels)
                    chanlabels = {chanlabels};
                end
                
                
                for i = 1:nchans
                    
                    if length(chanlabels)==nchans
                        ichanlabel = chanlabels{i};
                    else
                        ichanlabel = chanlabels;
                    end
                    if length(chanvals)==nchans
                        ichanval = chanvals(i);
                    else
                        ichanval = i;
                    end
                    if length(samplerates)==nchans
                        isamplerate = samplerates(i);
                    else
                        isamplerate = samplerates;
                    end
                    if length(tstarts)==nchans
                        itstart = tstarts(i);
                    else
                        itstart = tstarts;
                    end
                    if length(tends)==nchans
                        itend = tends(i);
                    else
                        if strcmp(isamplerate,'event')
                            itend = max(data_in);
                        else
                            itend = itstart + (length(data_in(:,i))-1)/isamplerate;
                        end
                    end
                    if length(unitss)==nchans
                        iunits = unitss(i);
                    else
                        iunits = unitss;
                    end
                    
                    
                    if ~isempty(data_in)
                        obj(i).data = data_in(:,i);
                    end
                    obj(i).chanlabel = ichanlabel;
                    obj(i).chanval = ichanval;
                    obj(i).samplerate = isamplerate;
                    obj(i).tstart = itstart;
                    obj(i).tend = itend;
                    obj(i).units = iunits;
                end
            end
            
        end
        
        % Methods
        function chanind = datchanind(obj,chans)
            % Get the index for the specified channel numbers or labels
            if isnumeric(chans)
                for i = 1:length(chans)
                    if find([obj.chanval]==chans(i))
                        chanind(i) = find([obj.chanval]==chans(i));
                    end
                end
            elseif ischar(chans) % Chanlabel specified (one)
                chanind = find(strcmp({obj.chanlabel},chans));     
            elseif iscell(chans) % Chanlabel specified (many)
                for i = 1:length(chans)
                    if find(strcmp({obj.chanlabel},chans{i}))
                        chanind(i) = find(strcmp({obj.chanlabel},chans{i}));
                    end
                end
            end
            chanind = chanind(chanind>0);
        end
        function datout = datchan(obj, chans)
            % Get a dat structure with just the specificied channels
            chaninds = datchanind(obj, chans);
            datout = obj(chaninds);
        end
        
        function datout = datdelete(obj,chans)
            chaninds = datchanind(obj,chans);
            datout = obj;
            datout(chaninds) = [];
        end
        
        function data = datchandata(obj,chan)
            % data = datchandata(obj,chan)
            % Get just the raw data from one channel
            
            if ~exist('chan','var') || length(obj)==1
                data = obj.data;
            else
                
                chanind = datchanind(obj,chan);
                datout = obj(chanind);
                if ~isempty(datout)
                    data = datout.data;
                else
                    error('No data or channel not found')
%                     data = [];
                end
            end
        end
        function datout = datseg(obj, seg, npoints)
            % datout = datseg(obj,[tstart tstop])
            % Return dat object with only the given time segment of data
            
            if ~exist('npoints','var')
                npoints = [];
            end
            
            if size(seg,2) ==1
                seg = seg';
            end
            
            datout = obj;
            
            for i = 1:length(obj)
                
                datout(i).data = [];
                
                for j = 1:size(seg,1)
                    
                    if ~isnumeric(obj(i).samplerate)
                        mask = obj(i).data>seg(j,1) & obj(i).data<seg(j,2);
                        datout(i).data = obj(i).data(mask);
                        datout(i).tstart = seg(j,1);
                        datout(i).tend = seg(j,2);
                        if strcmp(obj(i).samplerate,'event')
                            continue;
                        else
                            datout(i).samplerate = datout(i).samplerate(mask);
                        end
                        
                    else
                        
                        if isempty(npoints)
                            n = round((diff(seg(j,:))+2*eps)*obj(i).samplerate);
                        else
                            n = npoints;
                        end
                        
                        %                     startindex = round(obj(i).tstart + obj(i).samplerate*seg(1) + 1);
                        startindex = round((seg(j,1)-obj(i).tstart)* obj(i).samplerate + 1);
                        
                        %                     t = dattime(obj(i));
                        %                     startindex = find(t>=seg(1)-1e-4,1);
                        % Base sampling on total duration - so segments with the same duration
                        % and sample rate have the same length output
                        
                        stopindex = startindex + n-1;
                        
                        tstart = seg(j,1);
                        tend = seg(j,2);
                        if startindex < 1 || stopindex > length(obj(i).data)
                            %                         error('Time range exceeds data limits')
                            
                            
%                             if length(obj)==1
                                if startindex < 1
                                    startindex = 1;
                                    tstart = obj(i).tstart;
                                end
                                if stopindex > length(obj(i).data)
                                    stopindex = length(obj(i).data);
                                    tend = obj(i).tend;
                                end
%                             else
%                                 return
%                             end
                        end
                        %                         datout(i).data = obj(i).data(startindex:stopindex);
                        %                         datout(i).tstart = tstart;%t(startindex);
                        %                         datout(i).tend = tend;%t(stopindex); 11/1/13
                        
                        % Allow concatenation of multiple segments
                        datout(i).data = [datout(i).data; obj(i).data(startindex:stopindex)];
                        if j==1
                            datout(i).tstart = tstart;
                            datout(i).tend = tend;
                        else
                            datout(i).tend = datout(i).tend  + tend-tstart;
                        end
                        
                        
                    end
                end
            end
        end
        function data = datsegdata(obj, tstart, duration, alignsegs)
            % data = datsegdata(obj, tstart, duration, (alignsegs))
            %
            % Return the raw segmented data in a mtrix form with one
            % row per sample and one column per segment
            % since all segs must be the same length (for cont data), give tstart and
            % duration as input args
            %
            % For event data, returns a vector of all the event times,
            % aligned to each tstart
            %
            % DURATION can have two elements - [tbefore and tafter]
            
            if ~exist('alignsegs','var')
                alignsegs = 1;
            end
            if length(obj)>1
                error('Only input dat object with one channel - use datChan(obj, chan)')
            end
            
            if length(duration) == 2
                tbefore = duration(1);
                tafter = duration(2);
                duration = tafter - tbefore;
            elseif length(duration)==1
                tbefore = 0;
                tafter = duration;
            else
                tbefore  = 0;
                tafter = duration;
            end
            
            npoints = round( (tafter-tbefore)*obj(1).samplerate);
            
            data = [];
            for i = 1:length(tstart)
                                    
                if isnan(tstart(i))
                    data = [data NaN(npoints,1)];

                    continue
                end
                if obj.tend < tstart(i)+tafter || obj.tstart>tstart(i)+tbefore
                %*** HP 1/10/14
%                     warning('Excluding segments that fall outside range')
%                     fprintf('trial = %g s, tbefore = %g, tafter = %g, obj.tstart = %g, obj.tend = %g\n\n',tstart(i), tbefore, tafter,obj.tstart, obj.tend)
                    continue
                end
                if length(tafter)>1
                datTemp = datseg(obj, [tstart(i)+tbefore tstart(i)+tafter(i)], npoints);
                else
                datTemp = datseg(obj, [tstart(i)+tbefore tstart(i)+tafter], npoints);
                end

                if ~isnumeric(obj.samplerate)
                    if alignsegs
                        data = [data; datTemp.data  - tstart(i)];
                    else
                        data = [data; datTemp.data];
                    end
                    else
                    data = [data datTemp.data(:)];
                end
                
            end
            
        end
        
        function [datout, ncycles] = datsegmeans(obj, tstart, segment, shiftsegs,nonan, binwidth)
            % datout = datsegmeans(obj, tstart, duration, shiftsegs,nonan)
            % Returns mean values for all segments in a dat structure
            % tstart - time of "triggers" to align to
            % segment - can be one or two elements, specificy the time before and
            % after tstart
            % shiftsegs = 1 to shift all data so first point is 0
            % nonan = 1 ignores all segments containing any NaNs
            
            if ~exist('shiftsegs','var')
                shiftsegs = 0;
            end            
            if ~exist('nonan','var')
                nonan = 0;
            end            
            if ~exist('binwidth','var')
                binwidth = 1/obj(1).samplerate;
            end
            if numel(segment)==2
                startTime = segment(1);
                duration = diff(segment);
            else
                startTime = 0;
                duration = segment;
            end
            
            for i = 1:length(obj)
                datout(i) = obj(i);

                if isempty(obj(i).data)
                    continue;
                end
                
                if ~isnumeric(obj(i).samplerate)
                    segData = [];
                                       
                    nbins = round(duration/binwidth);
                    for j = 1:length(tstart)
%                         bins = (0:nbins-1)*binwidth + tstart(j) + startTime;
                        bins = (0:nbins)*binwidth + tstart(j) + startTime;
                        counts = histc(obj(i).data,bins);
                        segData(:,j) = counts(1:end-1)/binwidth;                        
                    end
                    datout(i).samplerate = 1/binwidth;
                    
                else
                    segData = datsegdata(obj(i), tstart, segment);
                    if shiftsegs % set first data point to zero
                        offsets = repmat(segData(1,:),[size(segData,1), 1]);
                        segData = segData-offsets;
                    end
                    if nonan
                        nansegs = find(sum(isnan(segData),1));
                        segData(:,nansegs)=[];
                    end
                end                
                datout(i).data = nanmean(segData,2);
                ncycles = size(segData,2);                
                datout(i).tstart = startTime;
                datout(i).tend = startTime+duration;
             
            end            
        end
        
        function [datout, ncycles] = datsegsems(obj, tstart, segment, shiftsegs,nonan)
            % datout = datsegmeans(obj, tstart, duration, shiftsegs,nonan)
            % Returns mean values for all segments in a dat structure
            % Must specify tstart and duration
            % shiftsegs = 1 to shift all data so first point is 0
            % nonan = 1 ignores all segments containing any NaNs
            
            if ~exist('shiftsegs','var')
                shiftsegs = 0;
            end
            
            if ~exist('nonan','var')
                nonan = 0;
            end
            
              
            for i = 1:length(obj)
                datout(i) = obj(i);

                if isempty(obj(i).data) || ~isnumeric(obj(i).samplerate)
                    continue;
                end
                
                segData = datsegdata(obj(i), tstart, segment);
                if shiftsegs % set first data point to zero
                    offsets = repmat(segData(1,:),[size(segData,1), 1]);
                    segData = segData-offsets;
                end
                if nonan
                    nansegs = find(sum(isnan(segData),1));
                    segData(:,nansegs)=[];
                end
                
                datout(i).data = sem(segData,2);
                ncycles = size(segData,2);                
                if numel(segment)==1
                    datout(i).tstart = 0;
                    datout(i).tend = segment;
                else
                    datout(i).tstart = segment(1);
                    datout(i).tend = segment(2);
                end
            end            
        end
        
        
        
        function datout = dateventtocont(obj, npointsorsamplerate, tstart, tend)
            %    dateventtocont(obj, samplerate) Uses tstart and tend from
            %    obj
            %    dateventtocont(obj, npoints, tstart, tend)
            % Make sure units of the event channel are in seconds!
            
            events = obj.data;
           
            if nargin<3
                fs = npointsorsamplerate;
                tstart = obj.tstart;
                tend = obj.tend;
                npoints = round((obj.tend-obj.tstart)*fs+1);
            elseif nargin<4
                fs = npointsorsamplerate;
                tend = obj.tend;
                npoints = round((tend-tstart)*fs+1);

            else
                npoints = npointsorsamplerate;
                fs = (npoints-1)/(tend-tstart);
            end
            if any(events<tstart); events = events(events>tstart); warning('dateventtocont events before tstart time dropped'); end
            currdata = zeros(1,npoints);                        
            currdata(round((events-tstart)*fs+1))=1;            
                        
            datout = dat(currdata, obj.chanlabel, obj.chanval, fs, tstart, tend, obj.units);
            
        end
        function datrange = datrange(obj,chan)
            % datRange Returns range (y-span) of data in specified channel
            %   chan is optional
            
            if exist('chan','var')
                datCurr = datchan(obj, chan);
                if numel(datCurr)~=1
                    error('Only pick one channel.')
                end
            elseif numel(obj)==1
                datCurr = obj;
            else
                error('Only pick one channel.')
            end
            datrange = [min(datCurr.data) max(datCurr.data)];
        end
        function datout = smooth(obj,varargin)
            % Smooths data in the dat structure DATIN
            % DATOUT = datsmooth(DATIN, [SPAN], [CHANS])
            %   DATIN - required input dat structure
            %   SPAN  - optional number of points to smooth
            %   CHANS - optional subset of channels to smooth (values
            %   or cell array of labels)
            
            
            p = inputParser;
            
            addRequired(p,'obj',@(x)isa(x,'dat'));
            addOptional(p,'span',5,@isnumeric);
            addOptional(p,'chans',[obj.chanval],@(x) isnumeric(x) || ischar(x) || iscell(x));
            
            parse(p,obj,varargin{:})
            
            obj = p.Results.obj;
            span = p.Results.span;
            chans = p.Results.chans;
            
            chanInd = datchanind(obj,chans);
            
            datout = obj;
            for i = 1:length(chanInd)
                if isnumeric(obj(chanInd(i)).samplerate)
                    datout(chanInd(i)).data = smooth(obj(chanInd(i)).data, span,'moving');
                end
            end
        end
        function datout = datsmooth(obj,varargin)
            % Smooths data in the dat structure DATIN
            % DATOUT = datsmooth(DATIN, [SPAN], [CHANS])
            %   DATIN - required input dat structure
            %   SPAN  - optional number of points to smooth
            %   CHANS - optional subset of channels to smooth (values
            %   or cell array of labels)
            
           datout = smooth(obj,varargin{:});
        end
        function datout = datbin(obj, binwidth)
            % Returns a dat structure with binned data
            % Specify binwidth in seconds
            
            databinned = zeros(1,floor(length(obj.data)/(binwidth*obj.samplerate)));
            
            for i = 1:length(databinned)
                istart = round(1 + (i-1)*binwidth*obj.samplerate);
                istop = round(i*binwidth*obj.samplerate);
                databinned(i) = nanmean(obj.data(istart:istop));
                
            end
            
            datout = obj;
            datout.data = databinned;
            datout.samplerate = 1/binwidth;
        end
        function datout = downsample(obj, varargin)
            % datDownsample Downsample input signal
            %   Downsample by keeping every nth sample starting with the
            %   first
            
            p = inputParser;
            
            addRequired(p,'obj',@(x)isa(x,'dat'));
            addOptional(p,'n',10,@isnumeric);
            
            parse(p,obj,varargin{:})
            
            obj = p.Results.obj;
            n = p.Results.n;
                             
            datout = obj;
            for i = 1:length(obj)
                if isnumeric(obj(i).samplerate)                                     
                        datout(i).data = downsample(obj(i).data,n);
                        datout(i).data = datout(i).data(:);
                        datout(i).samplerate = datout(i).samplerate/n;                                   
                end
            end
        end
        
         function datout = upsample(obj, n)
            % upsample Upsample input signal                           
            
            datout = obj;
            for i = 1:length(obj)
                if isnumeric(obj(i).samplerate)
                    
                    datout(i).data = resample(obj(i).data,n,1);
                    datout(i).data = datout(i).data(:);
                    datout(i).samplerate = datout(i).samplerate*n;
                end
            end
         end
        
         function datout = interp(obj, fs)
            %  interp Interpolate to new samplerate
            %  datout = interp(obj, fs)
            %  datout = interp(obj, ttNew)
                                   
            datout = obj;
            for i = 1:length(obj)
                if isnumeric(obj(i).samplerate)
                    ttOld = dattime(obj(i));
                    nPointsNew = (obj(i).tend - obj(i).tstart)*fs;
              
                    if length(fs)>1
                        ttNew = fs;
                    else
                        ttNew = obj(i).tstart + (0:(nPointsNew-1))/fs;
                    end
                    
                    datout(i).data = interp1(ttOld, obj(i).data,ttNew);                   
                    datout(i).data = datout(i).data(:);
                    datout(i).samplerate = 1/mean(diff(ttNew));
                    datout(i).tend = ttNew(end);
                end
            end
         end
        
         
         function datout = datlowpass(datin, fc)
             % datlowpass(datin, fc)
             % Filter dat channel with a lowpass butterworth filter
             % Use 100 Hz for eye vel filter
             
             datout = datin;
             for i= 1:length(datin)
                 datcur = datin(i);
                 
                 if strcmp(datcur.samplerate,'event')
                     continue;
                 end
                 
                 datcur.data = double(datcur.data);
                 fs = datcur.samplerate;
                 N = 3;      % Filter order
                 [b,a] = butter(N, fc/fs);
                 dataFiltered=filtfilt(b,a,double(datcur.data));
                 datout(i).data = dataFiltered;
             end
         end
         
        function datout = resettime(obj,time0)
            if ~exist('time0','var')
                time0=0;
            end
            datout = obj;
            for i = 1:length(obj)
                
                
                if ~isnumeric(obj(i).samplerate)
                    datout(i).data = obj(i).data - obj(i).tstart  + time0;
                end
                
                datout(i).tend = obj(i).tend - obj(i).tstart + time0;      
                datout(i).tstart = time0;
    
                
            end
        end
        
        function t = dattime(obj, seg)
            % t =  DATTIME(obj)
            % t = DATTIME(obj,seg)
            %
            % obj must be single channel
            
            npoints = length(obj(1).data);
            t = (0:npoints-1)/obj(1).samplerate + obj(1).tstart;            
            if exist('seg','var')
                mask = t>=seg(1)-eps & t<seg(2)-eps;
                t(~mask) = [];
            end
            t = t(:);
        end
        
        
        function data = dattimeind(obj,t)            
            % data = dattimeind(obj,t)
            % Returns data point at specified time. obj must be 1 channel
            datt = dattime(obj);
            for i = 1:length(t)
            [~,temp] = min(abs(t(i)-datt));
            ind(i) = temp(1);
            end
            data = obj.data(ind);
        end
        
        function [freq, pwr] = datfft(obj)
            
            n = length(obj.data);
            nhalf = floor(n/2);
            
            x = fft(obj.data);
            pwr = abs(x(1:nhalf));
            freq = (0:nhalf-1)'*obj.samplerate/n;
        end
        function varargout = plot(obj, varargin)
            %  datPlot   Simple plotting function.
            %   datPlot(obj) plots all channels in the dat structure
            %
            %   datPlot(obj, [tstart tend], 'Range',[ymin ymax],'Overlaid',1)
            
            p = inputParser;
            
            addRequired(p,'obj',@(x)isa(x,'dat'));
            addOptional(p,'seg',[]);
            addParamValue(p,'Range',[]);
            addParamValue(p,'Overlaid',0);
            addParamValue(p,'Color',[]);
            addParamValue(p,'Clipping','on');
            
            parse(p,obj,varargin{:});
            
            obj = p.Results.obj;
            seg = p.Results.seg;
            plotrange = p.Results.Range;
            overlaid = p.Results.Overlaid;
            color = p.Results.Color;
            clipping = p.Results.Clipping;
            
            emptychans = cellfun(@isempty,{obj.data});
            
            if all(emptychans)
                return
            end
            
            obj = obj(~emptychans);
            
            if isempty(seg)
                datdisp = obj;
            else
                datdisp = datseg(obj, seg);
            end
            
            
            
            if numel(plotrange)==1
                plotrange = [-plotrange, plotrange];
            end
            
            nchans = length(datdisp);
            if isempty(color)
                colors = hsv(nchans);
            elseif length(color)==1
                colors = repmat(colorspec(color),nchans,1);
            elseif ~iscell(color)
                colors = color;
                if size(color,1)==1
                    colors = repmat(color,nchans,1);
                end
                
            else
                for i = 1:length(color)
                    colors(i,:) = colorspec(color{i});
                end
            end
            lightmask = sum(colors,2)>1.5;
            colors(lightmask,:) = max(0,colors(lightmask,:)-.5);   % Darken
            
            h = [];
            %% Main plotting loop
            for i = 1:nchans
                
                if isempty(datdisp(i))
                    continue
                end
                
                color = colors(i,:);
                
                % Invert the channel number for plotting
                inv = nchans-i+1;
                
                % Set up axis labels and limits
                if ~overlaid
                    haxis(inv) = subplot(nchans,1,i);
                    
                    s1 = get(gca,'position');
                    ylabel(datdisp(inv).chanlabel,'FontWeight','normal','FontSize',10);
                    %ylabel([datdisp(inv).chanlabel ' (' datdisp(inv).units ')'])
                    %                     set(get(gca,'YLabel'),'FontWeight','normal','FontSize',10)
                    set(gca,'FontWeight','normal','FontSize',10)
                    xlim([datdisp(inv).tstart datdisp(inv).tend])
                    if  ~isnumeric(datdisp(inv).samplerate)
                        ylim([-1 1]);
                    elseif ~isempty(plotrange)
                        ylim(plotrange);
                    else
                        if islogical(datdisp(inv).data)
                            ylim([-.1 1.1])
                        elseif all(ismember(datdisp(inv).data,[-1 0 1]))
                            ylim([-1.1 1.1])
                        else
                            ylims = [prctile(datdisp(inv).data,.5) prctile(datdisp(inv).data,99.5)];
                            if ylims(2) > ylims(1)
                                ylim(ylims)
                            end
                        end
                    end
                    box off
                    set(gca,'position',s1);
                    if i<nchans
                        set(gca,'XColor','w')
                    else
                        xlabel('Time (s)')
                    end
                    
                else
                    xlabel('Time (s)')
                    haxis(1) = gca;
                    hold on;
                end
                
                
                % For normal time serieschannels
                if isnumeric(datdisp(inv).samplerate)
                    t = dattime(datdisp(inv));
                    dispdata = datdisp(inv).data;
                    
                    h(inv) = line(t,dispdata,'LineWidth',1,'Color',color);
                    
                else % For event channels
                    y = zeros(1,length(datdisp(inv).data));
                    if ~isempty(y)
                        h(inv) = line(datdisp(inv).data,y,'Marker','+','LineStyle','none','Color',color);
                    end
                end
                
            end
            
            %% Allow interactive zooming on all subplots
            dragzoom
            linkaxes(haxis,'x');
            
            set(h,'Clipping',clipping)
            
            %% Define output handle
            if nargout==1
                varargout = {h};
            end
            
        end % function datplot
        function datplot(obj,varargin)
            plot(obj,varargin);
        end
        
        
        function datout = deriv(obj)
            datout = obj;
            datout.data = diff(obj.data)*obj.samplerate;
            datout.units = [obj.units '/s'];
        end
        
        % Standard arithmetic functions
        function num = mean(obj)
            for i = 1:length(obj)
                num(i) = mean(obj(i).data);
            end
        end
        
        function datout = datmean(obj)
            data = NaN(length(obj(1).data),length(obj));
            for i = 1:length(obj)
                if ~isempty(obj(i).data)
                data(:,i) = obj(i).data;
                end
            end
            datout = obj(1);
            datout.data = nanmean(data,2);
        end
        
        function datout = plus(obj,a)
            datout = obj;
            for i = 1:numel(obj)
                if isa(a,'dat')
                    if length(a) > 1
                        num = a(i).data;
                    else
                        num = a.data;
                    end
                else
                    num = a;
                end
                datout(i).data = obj(i).data + num;
            end
        end
        function datout = minus(obj,a)
            datout = plus(obj,-a);
        end
        function datout = uminus(obj)
            datout = obj;
            for i = 1:numel(obj)
                if isnumeric(obj(i).samplerate);
                    datout(i).data = -obj(i).data;
                end
            end
        end
        function datout = times(obj,a)
            if ~isa(obj,'dat')
                temp = a;
                a = obj;
                obj = temp;
            end
            datout = obj;
            for i = 1:numel(obj)
                if isa(a,'dat')
                    datout(i).data = obj(i).data.*a(i).data;
                else
                    datout(i).data = obj(i).data .* a;
                end
            end
        end
        function datout = mtimes(obj, a)
            datout = times(obj, a);
        end
        function datout = mrdivide(obj,a)
            datout = obj;
            for i = 1:numel(obj)
                if isa(a,'dat')
                    datout(i).data = obj(i).data ./ a(i).data;
                else
                    datout(i).data = obj(i).data ./ a;
                end
            end
        end
        function datout = mldivide(obj,a)
            datout = mrdivide(obj,a);
        end
        function datout = ldivide(obj,a)
            datout = mrdivide(obj,a);
        end
        function datout = rdivide(obj,a)
            datout = mrdivide(obj,a);
        end
        
        
    end % METHODS
end % CLASSDEF DAT
