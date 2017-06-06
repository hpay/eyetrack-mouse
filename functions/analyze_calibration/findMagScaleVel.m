function  [scaleCh1, scaleCh2] = findMagScaleVel(freq)
% Simple velocity calibration
% Fit sine wave to both video and magnet velocity data.
% Calibration factor is equal to ratio of Avid/Amag
% Channel with best r2 to sine wave fit is selected, other channel is set
% to 0
% Output: scaleCh1, scaleCh2 saved
%
% Updated 2/23/15 
% Hannah Payne <hpayne@stanford.edu>
% Raymond Lab

if ~exist('freq','var')
    freq = 1;
end

close all; 
magthresh = 10;
vidthresh = 200;

% SET ANGLE HERE
theta = 40;

%% Load magnet data from Spike2
pathname = cd;
fList = dir;

% Find file with "calib" and ".smr" somewhere in the name
temp = arrayfun(@(x) ~isempty(regexpi(x.name, '^.*calib.*\.smr$')),fList);
filenameroot = fList(temp).name;
filenameroot = filenameroot(1:end-4);
if ~isempty(regexpi(filenameroot ,'^.*calibR.*$')) %% Right eye
    channels = [4 7 8 10];%*** Right eye only  
    disp('Loading RIGHT EYE channels')
else
    channels = [4 5 6 10]; %** left eye
end

% Load spike2
magnet = importSpike(fullfile(pathname,[filenameroot '.smr']),channels);

%% Load video
A=load(fullfile(pathname, 'vid_results.mat'));
results = A.results;

[vidH, vidV, ~] = calcEyeAngle(results, theta);

tvid = (results.time1+results.time2)/2;
tvid = tvid-tvid(1);


% dt = mean(diff(tvid));
% vidH_vel = diff(vidH)/dt;
% [amp_vid_predesaccade,phase_vid_predesaccade] = fitsine(dt, vidH_vel(1:400), 1,1)

% % Alternate method using just one camera
% % rP = calcEyeRadius(results,theta);
% posH1 = real(asind((results.pupil1(:,1) - results.cr1a(:,1))./rP)); % in degrees cam 1 estimate
% posH2 = real(asind((results.pupil2(:,1) - results.cr2b(:,1))./rP)); % in degrees cam 2 estimate

%% Subselect and desaccade magnet
lightpulse = magnet(end).data;

if exist(fullfile(pathname,[filenameroot '.xlsx']),'file')
    seg(1) = xlsread(fullfile(pathname,[filenameroot '.xlsx']),'D2:D2');
else
    seg(1) = lightpulse(1);
end

seg(2) = seg(1)+tvid(end)

magnetSeg = resettime(datseg(magnet,seg));
samplerateM = magnetSeg(1).samplerate;
magnetSeg(1).data = smooth([diff(smooth(magnetSeg(1).data,25)); 0]*samplerateM,25); % Convert to velocity
magnetSeg(1).units = 'deg/s'; magnetSeg(1).chanlabel = 'hhvel';

% filter with 50 Hz low pass filter
mag1 = datlowpass(magnetSeg(2),50); 
mag2 = datlowpass(magnetSeg(3),50);

tmag = dattime(mag1);

mag1Vel = [diff(mag1.data); 0]*samplerateM;
mag1_saccademask = ~isnan(desaccade(mag1Vel,samplerateM,magthresh,0)); 
mag1Velplot = smooth(mag1Vel,50);  
h = 1:3;
figure(h(1)); clf; title('Magnet 1 desaccading'); hold on;  

plot(tmag, mag1Velplot,'Color',.8*[1 1 1]); mag1Velplot(~mag1_saccademask) = NaN; 
plot(tmag, mag1Velplot,'r'); ylim([-5 5]); 

mag2Vel = [diff(mag2.data); 0]*samplerateM;
mag2_saccademask = ~isnan(desaccade(mag2Vel,samplerateM,magthresh,0));
mag2Velplot = smooth(mag2Vel,50); 
figure(h(2)); clf; title('Magnet 2 desaccading'); hold on; 
plot(tmag, mag2Velplot,'Color',.8*[1 1 1]);  mag2Velplot(~mag2_saccademask) = NaN;
plot(tmag, mag2Velplot,'r');  ylim([-5 5]); 


%% Desaccade video
vidH_upsample = interp1(tvid,vidH,tmag(:),'linear');

% Keep the big NaNs NaN
gapCumsum = cumsum(isnan(vidH_upsample));
gapDiff = [0; diff(isnan(vidH_upsample))];
gapIndex = cumsum(gapDiff).*cumsum(gapDiff>0);
maxGap = 250; % 250 ms
gapLength = diff([0; gapCumsum(gapDiff==-1)]);

% Fill in short NaNs
vidH_upsample = inpaint_nans(vidH_upsample);

% Keep the big NaNs
for itemp = 1:length(gapLength)
    if gapLength(itemp) > maxGap
    vidH_upsample(gapIndex==itemp) = NaN;
    end
end


vidVel = [diff(vidH_upsample); 0]*samplerateM;
vid_saccademask = ~isnan(desaccade(vidVel,samplerateM,vidthresh,0)); 

tVidGood = tvid(~isnan(vidH));
vid_saccademask(tmag>tVidGood(end)) = 0;

vidVelplot = smooth(inpaint_nans(vidVel),50,'moving'); 

figure(h(3)); clf; title('Video desaccading'); hold on
plot(tmag, vidVelplot,'Color',.8*[1 1 1]);  vidVelplot(~vid_saccademask) = NaN;
plot(tmag, vidVelplot,'r');  ylim([-50 50]); 


%% Average position of each
meanMag1Pos = nanmean(mag1.data(vid_saccademask & mag1_saccademask))
meanMag2Pos = nanmean(mag2.data(vid_saccademask & mag2_saccademask))
meanVidPos = nanmean(vidH_upsample(vid_saccademask))

%% Do sine fit on magnet and video sine waves
y1 = sin(2*pi*freq*tmag(:));
y2 = cos(2*pi*freq*tmag(:));
const = ones(size(y1));
vars = [y1 y2 const];

% ------------ Chair ------------
[bHead,bint,r,rint,stat] = regress(magnetSeg(1).data, vars);
headAmp = sqrt(bHead(1)^2+bHead(2)^2);
headPhase = rad2deg(atan2(bHead(2),bHead(1)));

% ------------ MAGNET 1------------
[bMag1,bint,r,rint,stat] = regress(mag1Vel(mag1_saccademask), vars(mag1_saccademask,:));
mag1Amp = sqrt(bMag1(1)^2+bMag1(2)^2);
mag1Phase = mod((rad2deg(atan2(bMag1(2),bMag1(1))) - headPhase),360)-180;
r2mag1 = stat(1);

% ------------ MAGNET 2------------
[bMag2,bint,r,rint,stat] = regress(mag2Vel(mag2_saccademask), vars(mag2_saccademask,:));
mag2Amp = sqrt(bMag2(1)^2+bMag2(2)^2);
mag2Phase = mod((rad2deg(atan2(bMag2(2),bMag2(1))) - headPhase),360)-180;
r2mag2 = stat(1);

% ------------ VIDEO ------------
[bVid,bint,r,rint,stat] = regress(vidVel(vid_saccademask), vars(vid_saccademask,:));
vidAmp = sqrt(bVid(1)^2+bVid(2)^2);
vidPhase = rad2deg(atan2(bVid(2), bVid(1)));
r2vid = stat(1);

%% save scale factor
scaleCh1_raw = vidAmp/mag1Amp * (2*(abs(mag1Phase)<90)-1); % Change sign if needed
scaleCh2_raw = vidAmp/mag2Amp * (2*(abs(mag2Phase)<90)-1); % Change sign if needed

if r2mag1>r2mag2
    scaleCh1 = scaleCh1_raw; 
    scaleCh2 = 0;
else
    scaleCh1 = 0;
    scaleCh2 = scaleCh2_raw; 
end
fprintf('scaleCh1 = %.2f, scaleCh2 = %.2f\n',scaleCh1, scaleCh2)
fprintf('r2   Ch1 = %.4f, r2   Ch2 = %.4f\n',r2mag1, r2mag2)
save(fullfile(cd, [filenameroot '.mat']),'scaleCh1', 'scaleCh2',...
    'vidAmp','mag1Amp','mag1Phase','mag2Amp','mag2Phase',...
    'r2mag1','r2mag2','r2vid','vidthresh','magthresh','freq',...
    'meanMag1Pos','meanMag2Pos','meanVidPos');

%% Plot fits
xlims = [0 min(50,tmag(end))];
figure(h(1)); hold on; plot(tmag,vars*bMag1,'k-'); box off;xlim(xlims)
a = get(gcf,'Position'); % Shrink vertically
set(gcf,'Position',[a(1) a(2) a(3) a(4)/2])
text(1, 4, sprintf('amp = %.3f, r2 = %.3f, scaleCh1 = %.2f',mag1Amp, r2mag1, scaleCh1_raw),'BackgroundColor','w')
saveas(gcf, 'fitMagnet1.png');

figure(h(2)); hold on; plot(tmag,vars*bMag2,'k-'); box off;xlim(xlims)
a = get(gcf,'Position'); % Shrink vertically
set(gcf,'Position',[a(1) a(2) a(3) a(4)/2])
text(1, 4, sprintf('amp = %.3f, r2 = %.3f, scaleCh2 = %.2f',mag2Amp, r2mag2, scaleCh2_raw),'BackgroundColor','w')
saveas(gcf, 'fitMagnet2.png');

figure(h(3)); hold on; plot(tmag, vars*bVid,'k-'); box off;xlim(xlims)
a = get(gcf,'Position'); % Shrink vertically
set(gcf,'Position',[a(1) a(2) a(3) a(4)/2])
text(1, 40, sprintf('amp = %.3f, r2 = %.3f',vidAmp, r2vid),'BackgroundColor','w')
saveas(gcf, 'fitVideo.png');


end