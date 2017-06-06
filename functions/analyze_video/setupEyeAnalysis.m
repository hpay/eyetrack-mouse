function setupEyeAnalysis
% setupEyeAnalysis
% Set the settings for dual-angle video-oculography
%
% 2014 
% Hannah Payne <hpayne@stanford.edu>
% Raymond Lab, Stanford University 

%% DEFAULT PARAMETERS
% Select initial pupil location & radius automatically or manually
p.manual = 0;       % {0,1}

% Pupil radius initial parameters
p.radiiPupil = 30;  % (px) default radius 

% corneal refleaction radii bounds
p.radiiCR = [4 8];  % (px) lower and upper limits for CR radius
p.CRthresh = 10;  % (default 10) Threshold for subtracting background (default 10)   
p.CRfilter = 3; % (Optional, default is 8, minimum is 3)
%               The radius of the filter used in the search of local
%               maxima in the accumulation array. To detect circles whose
%               shapes are less perfect, the radius of the filter needs
%               to be set larger.

% Fraction features needed for pupil detection - smaller = faster
p.minfeatures = .6; % (fraction, [0 1])

% Improve image contrast 
p.imAdjust = 1;     % 0 is faster, 1 may work better

if exist(fullfile(cd,'vid_settings.mat'),'file');
    load(fullfile(cd,'vid_settings.mat'));
    p.radiiPupil = round(mean(p.radiiPupil));
end

%% ALLOW PARAMETER ADJUSTMENT AND TEST
ok = 'No';
while strcmp(ok,'No')   
    % Get input
    prompt = {'manual:', 'pupil radius guess:','CR radius guess:','minfeatures (0-1):','imAdjust:'};
    dlg_title = 'Settings for eyeAnalysis';
    defAns =   cellfun(@num2str,{p.manual, p.radiiPupil, p.radiiCR, p.minfeatures, p.imAdjust},'UniformOutput',0);
    answerStr = inputdlg(prompt,dlg_title,1,defAns);
    answerNum = cellfun(@str2num, answerStr,'UniformOutput',0);
    [p.manual, p.radiiPupil, p.radiiCR, p.minfeatures, p.imAdjust]= answerNum{:};
    
    %% Run test of setttings and save them
    videoTest(p)
    
    %% Ask if ok
    ok = questdlg('Accept settings?','Accept settings?','Yes','No','Yes');
end



