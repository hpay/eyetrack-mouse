function d = readSpikeFile(filename,chans)

%% readSpikeFile.m imports Cambridge Electronic Design Spike2 files
%
% Requires the MATLAB SON library version 2.4 or higher.
%
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 07/06
% Copyright © The Author & King's College London 2006-2007
% -------------------------------------------------------------------------
%
% Revisions:
% Modified by Hannah Payne 07/2012 for Raymond lab
% If no channels specified, output is structure c with channel list


if SONVersion('nodisplay')<2.31
    errordlg('ImportSMR: An old version of the SON library is on the MATLAB path.\nDelete this and use the version in sigTOOL');
    which('SONVersion');
    return;
end

[pathname, filename2, extension]=fileparts(filename);

fprintf('Reading Spike2 file: %s\n',fullfile(pathname, filename2));

if strcmpi(extension,'.smr')==1 || strcmpi(extension,'.srf')==1
    % Spike2 for Windows source file so little-endian
    fid=fopen(filename,'r','l');
elseif strcmpi(extension,'.son')==1
    % Spike2 for Mac file
    fid=fopen(filename,'r','b');
else
    warning('%s is not a Spike2 file\n', filename);
    return
end

if fid<0
    errordlg('File not found check path')
    return
end


% get list of valid channels
F=SONFileHeader(fid);
c=SONChanList(fid);

%% Display all the channels
for i = 1:length(c)
    fprintf('%i %s\n',c(i).number,c(i).title);
end
disp(' ')

%% Import the data.
if ~exist('chans','var') || isempty(chans)
    d = c;
else
    
    index = 1;
    for i=1:length(c)
        
        
        chan=c(i).number;
        
        if ismember(chan, chans)
            
            % For each channel, call the SON library function SONGetChannel
            try
                [data,header]=SONGetChannel(fid, chan);
                
                if isempty(data)
                    % Empty channel
                    warning(sprintf('Empty channel %i',chan))
                    continue
                end
                
            catch msgid
                % Failed, go to next channel
                warning('Channel %g failed.\n',chan);
                continue;
            end
            
            
            hdr.channel=chan;
            hdr.source=dir(header.FileName);
            hdr.source.name=header.FileName;
            hdr.title=header.title;
            hdr.comment=header.comment;
            if strcmpi(hdr.title,'Keyboard')
                hdr.markerclass='char';
            else
                hdr.markerclass='uint8';
            end
            
            switch header.kind
                case {1,9}% Waveform int16 or single in SMR file
                    
                    imp.tim(:,1)=header.start;
                    imp.tim(:,2)=header.stop;
                    imp.adc=data;
                    imp.mrk=zeros(size(imp.tim,1),4,'uint8');
                    
                    if size(imp.adc,2)==1
                        hdr.channeltype='Continuous Waveform';
                        hdr.channeltypeFcn='';
                        hdr.adc.Labels={'Time'};
                    else
                        hdr.channeltype='Episodic Waveform';
                        hdr.adc.Labels={'Time' 'Epoch'};
                    end
                    
                    hdr.adc.TargetClass='adcarray';
                    hdr.adc.SampleInterval=[header.sampleinterval 1e-6];
                    if header.kind==1
                        hdr.adc.Scale=header.scale/6553.6;
                        hdr.adc.DC=header.offset;
                    else
                        hdr.adc.Scale=1;
                        hdr.adc.DC=0;
                    end
                    hdr.adc.Func=[];
                    hdr.adc.Units=header.units;
                    hdr.adc.Multiplex=header.interleave;
                    hdr.adc.MultiInterval=[0 0];%not known from SMR format
                    hdr.adc.Npoints=header.npoints;
                    hdr.adc.YLim =[header.min header.max]*hdr.adc.Scale+hdr.adc.DC;
                    
                    hdr.tim.Class='tstamp';
                    hdr.tim.Scale=F.usPerTime;
                    hdr.tim.Shift=0;
                    hdr.tim.Func=[];
                    hdr.tim.Units=F.dTimeBase;
                    
                case {2,3}% Event+ or Event- in SMR file
                    imp.tim(:,1)=data;
                    imp.adc=[];
                    imp.mrk=zeros(size(imp.tim,1),4,'uint8');
                    if header.kind==2
                        hdr.channeltype='Falling Edge';
                    else
                        hdr.channeltype='Rising Edge';
                    end
                    hdr.channeltypeFcn='';
                    hdr.adc=[];
                    
                    hdr.tim.Class='tstamp';
                    hdr.tim.Scale=F.usPerTime;
                    hdr.tim.Shift=0;
                    hdr.tim.Func=[];
                    hdr.tim.Units=F.dTimeBase;
                    
                case {4}% EventBoth in SMR file
                    if header.initLow==0 % insert a rising edge...
                        data=vertcat(-1, data);   % ...if initial state is high
                    end
                    
                    % EDIT HP to avoid losing assymetric pulses
                    if mod(length(data),2) ~= 0
                        data(end+1) = NaN;                        
                    end
                    imp.tim(:,1)=data(1:2:end);% rising edges
                    %             imp.tim(:,1)=data(1:2:end-1);% rising edges
                    imp.tim(:,2)=data(2:2:end);% falling edges
                    imp.adc=[];
                    imp.mrk=zeros(size(imp.tim,1),4,'uint8');
                    
                    hdr.channeltype='Pulse';
                    hdr.channeltypeFcn='';
                    hdr.adc=[];
                    
                    hdr.tim.Class='tstamp';
                    hdr.tim.Scale=F.usPerTime;
                    hdr.tim.Shift=0;
                    hdr.tim.Func=[];
                    hdr.tim.Units=F.dTimeBase;
                    
                case {5}% Marker channel in SMR file
                    imp.tim(:,1)=data.timings;
                    imp.adc=[];
                    imp.mrk=data.markers;
                    
                    hdr.channeltype='Edge';
                    hdr.channeltypeFcn='';
                    hdr.adc=[];
                    
                    hdr.tim.Class='tstamp';
                    hdr.tim.Scale=F.usPerTime;
                    hdr.tim.Shift=0;
                    hdr.tim.Func=[];
                    hdr.tim.Units=F.dTimeBase;
                    
                case {6}% int16 ADC Marker in SMR file
                    imp.tim(:,1)=data.timings;
                    % 24.02.08 remove -1 and include interleave factor
                    imp.tim(:,2)=data.timings...
                        +(SONGetSampleTicks(fid,chan)*(header.preTrig));
                    imp.tim(:,3)=data.timings...
                        +(SONGetSampleTicks(fid,chan)*(header.values/header.interleave-1));
                    
                    imp.adc=data.adc;
                    imp.mrk=data.markers;
                    
                    hdr.channeltype='Framed Waveform (Spike)';
                    hdr.channeltypeFcn='';
                    
                    hdr.adc.Labels={'Time' 'Spike'};
                    hdr.adc.TargetClass='adcarray';
                    hdr.adc.SampleInterval=[header.sampleinterval 1e-6];
                    hdr.adc.Scale=header.scale/6553.6;
                    hdr.adc.DC=header.offset;
                    hdr.adc.YLim=[double(min(data.adc(:)))*hdr.adc.Scale+hdr.adc.DC...
                        double(max(data.adc(:)))*hdr.adc.Scale+hdr.adc.DC];
                    hdr.adc.Func=[];
                    hdr.adc.Units=header.units;
                    hdr.adc.Npoints(1:size(imp.adc,2))=header.values;
                    hdr.adc.Multiplex=header.interleave;
                    hdr.adc.MultiInterval=[0 0];%not known from SMR format
                    
                    hdr.tim.TargetClass='tstamp';
                    hdr.tim.Scale=F.usPerTime;
                    hdr.tim.Shift=0;
                    hdr.tim.Func=[];
                    hdr.tim.Units=F.dTimeBase;
                    
                case {7,8}% Real marker or text marker in SMR file
                    imp.tim(:,1)=data.timings;
                    switch header.kind
                        case 7
                            imp.adc=data.real;
                            hdr.channeltype='Edge';
                            hdr.adc.TargetClass='single';
                            hdr.channeltypeFcn='';
                            hdr.adc.Labels={'Single'};
                        case 8
                            imp.adc=data.text;
                            hdr.channeltype='Edge';
                            hdr.adc.TargetClass='char';
                            hdr.channeltypeFcn='SONMarkerDisplay';
                            hdr.adc.Labels={'Text'};
                    end
                    imp.mrk=data.markers;
                    hdr.adc.SampleInterval=NaN;
                    hdr.adc.Func=[];
                    hdr.adc.Scale=1;
                    hdr.adc.DC=0;
                    hdr.adc.Units='';
                    hdr.adc.Multiplex=NaN;
                    hdr.adc.MultiInterval=[0 0];%not known from SMR format
                    hdr.adc.Npoints(1:size(imp.adc,2))=header.values;
                    
                    hdr.tim.TargetClass='tstamp';
                    hdr.tim.Scale=F.usPerTime;
                    hdr.tim.Shift=0;
                    hdr.tim.Func=[];
                    hdr.tim.Units=F.dTimeBase;
                    
                otherwise
                    continue
            end
            
            d(index).data = imp;
            d(index).header = hdr;
            
            fprintf('Channel %g read.\n', chan)
            
            index = index + 1;
            clear('imp','hdr','data','header');
            
        end
        
    end
end
disp('Spike file read!')
fclose(fid);

