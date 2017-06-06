%  results = eyeAnalysis
% Run this to eye track on saved image files
% Adds to a results structure that already has time info
% to do - have ROI picked automatically based on pupil radius size, and
% move with pupil

%% Parameters - try changing if needed
clear; close all; clc
load vid_settings

% Shouldn't need to change unless setup changes drastically
edge_thresh1 = 45; % Initial gradient threshold of pupil edge detection for cam 1
edge_thresh2 = 45; % Initial gradient threshold of pupil edge detection for cam 2
plot_all = 1;
debug_on = 0;
plot_handles1 = [];
plot_handles2 = [];
plot_handles_result = [];
goal_perc_complete = 0;

%% Load results time file
temp = load(fullfile(cd,'vid_time.mat'));
tt = mean([temp.results.time1(:) temp.results.time2(:)],2); 
tt = tt - tt(1);
nt = length(tt);

v1 = VideoReader('img1.avi');
v2 = VideoReader('img2.avi');
n_images = v1.Duration * v1.FrameRate + 1;

results.pupil1 = NaN(nt,5);
results.pupil2 = NaN(nt,5);
results.cr1a = NaN(nt,3);
results.cr2a = NaN(nt,3);
results.cr1b = NaN(nt,3);
results.cr2b = NaN(nt,3);
results.sync1 = NaN(nt,1);
results.sync2 = NaN(nt,1);
results.time1 = tt;
results.time2 = tt;

tic
figure(2); clf; %set(1,'WindowStyle','docked')
figure(1); clf; %set(1,'WindowStyle','docked')


%% Start looping
for i = 1:nt-1
    
    % Check end of file
    if ~hasFrame(v1) || ~hasFrame(v2); break;   end
    
    % Load images
    [img1, img2] = readFramesGray(v1, v2);
    
    % Use the previous pupil location as a starting point
    if i~=1
        p.pupilStart1 = results.pupil1(i-1,:);
        p.pupilStart2 = results.pupil2(i-1,:);
    end
    
    try
        %% CAMERA 1
        plot_all = 1;
        if isempty(plot_handles1); figure(1); subplot(1,2,1); end
        [results.pupil1(i,:), results.cr1a(i,:), results.cr1b(i,:), points1, edge_thresh1, results.sync1(i), plot_handles1] = ...
            detectPupilCR(img1,edge_thresh1+4,p.pupilStart1, p, 'PlotOn',plot_all,'debugOn',debug_on,'PlotHandles', plot_handles1);
    catch msgid
        fprintf('Error in cam 1 img %i\n',i)
        msgid
        edge_thresh1 = 45;
    end
    
    try
        %% CAMERA 2
        plot_all = 1;
        if isempty(plot_handles2); figure(1); subplot(1,2,2); end
        [results.pupil2(i,:), results.cr2a(i,:), results.cr2b(i,:),points2, edge_thresh2, results.sync2(i), plot_handles2] = ...
            detectPupilCR(img2,edge_thresh2+4,p.pupilStart2, p, 'PlotOn',plot_all,'debugOn',debug_on, 'PlotHandles',plot_handles2);
        
    catch msgid
        fprintf('Error in cam 2 img %i\n',i)
        msgid
        edge_thresh2 = 45;
    end;
    drawnow;
    
    % Report % complete
    perc_complete = i/nt*100;
    if perc_complete >= goal_perc_complete
        
        fprintf('%1.f%% \n', perc_complete);
        goal_perc_complete = perc_complete + 1;
        
        % Save incrementally
        save('videoresults.mat','results')
        
        % Plot current results
        set(0,'CurrentFigure',2)
        plot_handles_result = plotResults(results, plot_handles_result); %set(gcf,'WindowStyle','docked');
    end
    drawnow
end

%% Plot and save
print(1, 'eyeAnalysis1.png')
print(2, 'eyeAnalysis2.png')

save('vid_results.mat','results')
fprintf('\nResults saved\n')
totalTime = toc;
fprintf('Processing time %f per image\n',totalTime/nt)

