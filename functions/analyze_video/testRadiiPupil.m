% testradiipupil
% Test pupil radius parameters for initial search
function [pupilStart1, val] = testRadiiPupil(img1, radiiPupil,dbgon)

%% Reserved parameters
alphaFRST = 1;           % Sensitivity to radial symmetry % 0.5
fracRemove = .15;

%% Camera
% Take the radial transform
imgRadialPupil = radialSymTransform(img1, radiiPupil, alphaFRST);
imgRadialPupil = removeBorder(imgRadialPupil, fracRemove);

% Find the lowest point (for dark pupil)
[pupilY1, pupilX1] = find(min(imgRadialPupil(:))==imgRadialPupil,1);
pupilStart1 = [pupilX1 pupilY1];

if dbgon
    figure; subplot(1,2,1); imagesc(imgRadialPupil); colormap(gray);hold on; plot(pupilStart1(1), pupilStart1(2),'+r'); axis image;
    subplot(1,2,2); imagesc(img1); colormap(gray);hold on; plot(pupilStart1(1), pupilStart1(2),'+r'); axis image;
    axis image;box off
end

val = min(imgRadialPupil(:)); % smaller is better for dark spots

end

