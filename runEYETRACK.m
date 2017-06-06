%% runEYETRACK
% Run a mouse video eye tracking calibration session
%
% NOTES
% Runs with two Logitech C310 cameras - test with supplied software first
%
% If you have already collected data or are using the example data,
% skip straight to analysis sections below (#3 and on)
%
% #1. View preview
% #2. Collect data
% #3. Set settings
% #4. Run video analysis
% #5. Run magnetic sensor calibration
%
% Dependencies: Image processing toolbox (to acquire images)
%
% 2013
% Hannah Payne
% Raymond Lab, Stanford University 
% hpayne@stanford.edu

% Add this folder and helper functions to path
addpath(fileparts(mfilename('fullpath')))
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'functions')))

options = {'Run calibration session','Analyze video data','Get calibration factor'};
option = questdlg('Choose an option:','runEyeTrack',options{1}, options{2}, options{3},options{1});

% Prompt to select path -- each calib. session should have its own folder
close all;
filepath = uigetdir;   cd(filepath)

switch option
        
        %% 1. Set up video cameras 
    case options{1}
        % Parameters for video recording (Logitech C310)
        P.fps = 30;           % (frames per second) approx. to estimate total time
        P.n_seconds = 120;    % (s) total time for calibration
        P.boxsize = 500;      % (px) size of initial ROI
                       
        % Set up videos
        close all; imaqreset; clc;
        [vid, P.res] = setupVideo(P.boxsize); pause(.2)
   
        % Show preview until figure closed
        [img1, img2] = eyecam(vid);
        
        % Select smaller ROI: 
        % 1. Click on camera image from the left side
        % 2. Click-drag to select rectangle
        vid = setROI(vid,img1, img2, P.res, P.boxsize);
        
        %% 2. Run the calibration session 
        % Start this at the same time (within 1 cycle) as the sine stim in Spike2
        close all
        eyecam(vid, 1, P.fps*P.n_seconds); % Will save img1, img2, and time.mat        
        
        %% 3. Set settings for video analysis
    case options{2}
        setupEyeAnalysis
        
        %% 4. Run video analysis
        eyeAnalysis     % locate pupil and CR
        
        % Save result figures
        print(1,'eyeAnalysis1.png','-dpng')
        print(2,'eyeAnalysis2.png','-dpng')       
        
        %% 5. Run calibration analysis
    case options{3}
        
        %% Based on velocity
        findMagScaleVel
        
end

