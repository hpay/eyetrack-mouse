function videoTest(p)
% Test the parameters for pupil and CR detetion, find correct radius
% automatically
%
% 2014 
% Hannah Payne <hpayne@stanford.edu>
% Raymond Lab, Stanford University 

% Paramters:
radiiPupils = 20:10:60; % (px) Range to test

%% Load image to test
close all

v1 = VideoReader('img1.avi');
v2 = VideoReader('img2.avi');

frame = 1;
img1 = readFrame(v1);
img2 = readFrame(v2);

img1 = rgb2gray(img1);
img2 = rgb2gray(img2);

%% Test pupil location
if length(p.radiiPupil)==1
    p.radiiPupil  = round([p.radiiPupil(1)*.8 p.radiiPupil(1)*1.2]); % Guess pupil radius in pixels
end

if p.manual
    [pupilStartLarge1, pupilStartLarge2] = testRadiiPupilManual(img1, img2, p);        
else
    
    % Set range of possible pupil radii hhere
    for i = 1:length(radiiPupils)
        radiiPupil = radiiPupils(i);
        [pupilStartLarge1all(:,i), val(1,i)] = testRadiiPupil(img1, round([radiiPupil*.8 radiiPupil*1.2]),1);
        [pupilStartLarge2all(:,i), val(2,i)] = testRadiiPupil(img2, round([radiiPupil*.8 radiiPupil*1.2]),1);
    end
    
    bestInd = find(mean(val)==min(mean(val)));
    pupilStartLarge1 = pupilStartLarge1all(:, bestInd);
    pupilStartLarge2 = pupilStartLarge2all(:, bestInd);
     
    close(setdiff(1:length(radiiPupils)*2, bestInd*2-1:bestInd*2))
    p.radiiPupil = radiiPupils(bestInd);
end

p.pupilStart1 = pupilStartLarge1(:)';
p.pupilStart2 = pupilStartLarge2(:)';

%% Start
edgeThresh1 = 35; % Initial gradient threshold of pupil edge detection for cam 1
edgeThresh2 = 35; % Initial gradient threshold of pupil edge detection for cam 2
plot_all = 1;
debug_on = 0;

ok = 0;
while ~ok && frame<300
    
    %% Select ROI  
    [img1, img2] = readFramesGray(v1, v2);
    frame  = frame+1;
        
    %% CAMERA 1
    try
        figure(1); subplot(1,2,1); cla;
        [pupil1, CR1a, CR1b] = detectPupilCR(img1,edgeThresh1,...
            p.pupilStart1, p,'PlotOn',plot_all,'DebugOn',debug_on); 
        ok = 1;
    catch msgid
        warning(msgid.message)
        ok = 0;
    end    
    
    %% CAMERA 2
    try
        figure(1); subplot(1,2,2); cla
        [pupil2, CR2a,CR2b] = detectPupilCR(img2,edgeThresh2, ...
            p.pupilStart2, p,'PlotOn',plot_all,'DebugOn',debug_on);
        continue; % If successful move on
    catch msgid
        warning(msgid.message)
        ok = 0;
    end
    
end


if p.manual
    mean_pupil_radius = mean([max(pupil1(3:4)) max(pupil2(3:4))])
else
    
    %% Run again using found pupil center and radius
    for k = 1:10
        
        mean_pupil_radius = mean([max(pupil1(3:4)) max(pupil2(3:4))])
        p.radiiPupil = round(mean_pupil_radius); % Guess pupil radius in pixels
        
        p.pupilStart1 = pupil1(1:2) ;
        p.pupilStart2 = pupil2(1:2) ;
        
        subplot(1,2,1); cla;
        [pupil1, CR1a, CR1b] = detectPupilCR(img1,edgeThresh1,p.pupilStart1, p,'PlotOn',plot_all,'DebugOn',debug_on); 

        subplot(1,2,2); cla;
        [pupil2, CR2a,CR2b] = detectPupilCR(img2,edgeThresh2, p.pupilStart2, p,'PlotOn',plot_all,'DebugOn',debug_on);

        drawnow;
        if sum(abs([pupil1(1:2) - p.pupilStart1  pupil2(1:2) - p.pupilStart2]))<10; break; end
        
    end
end

mean_CR_radius = nanmean([CR1a(3)  CR1b(3) CR2a(3) CR2b(3)])

%% Save settings
save('vid_settings','p')

end