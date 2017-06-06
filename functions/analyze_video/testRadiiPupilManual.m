% testradiipupil
% Test pupil radius parameters for initial search
function [pupilStart1,pupilStart2] = testRadiiPupilManual(img1, img2, p)

%% Reserved parameters
alphaFRST = 0.5;           % Sensitivity to radial symmetry

% %% Get images
% img1 = imread('img1.tiff','Index',1);
% img2 = imread('img2.tiff','Index',1);
    
if p.imAdjust
    img1 = imadjust(img1);
    img2  = imadjust(img2);
end

% Camera 1
clf; subplot(2,2,1); imagesc(img1);% colormap(gray);
axis image;box off

% Camera 2
subplot(2,2,2); imagesc(img2); %colormap(gray); 
axis image;box off



%% Get pupil center manually
disp('Click on each pupil center')
[x y] = ginput(2);
pupilX1 = x(1);
pupilX2 = x(2);
pupilY1 = y(1);
pupilY2 = y(2);


pupilStart1 = [pupilX1 pupilY1];
pupilStart2 = [pupilX2 pupilY2]

subplot(2,2,1)
hold on; plot(pupilStart1(1), pupilStart1(2),'+r')
subplot(2,2,2)
hold on; plot(pupilStart2(1), pupilStart2(2),'+r')




%% Select ROI
% pos1 = [pupilX1 pupilY1 0 0] + [-roiwidth/2 -roiheight/2 roiwidth roiheight];
% pos2 = [pupilX2 pupilY2 0 0] + [-roiwidth/2 -roiheight/2 roiwidth roiheight];
% pos1(1:2) = max(pos1(1:2),1);
% pos2(1:2) = max(pos2(1:2),1);
% pos1(3:4) = min(pos1(3:4),fliplr(size(img1))-pos1(1:2));
% pos2(3:4) = min(pos2(3:4),fliplr(size(img2))-pos2(1:2));
% pos1 = round(pos1);
% pos2 = round(pos2);
% img1 = img1(pos1(2):pos1(2)+pos1(4), pos1(1):pos1(1)+pos1(3));
% img2 = img2(pos2(2):pos2(2)+pos2(4), pos2(1):pos2(1)+pos2(3));


% pupilStart1 = round(pos1(3:4)/2);
% pupilStart2 = round(pos2(3:4)/2);
% 
% subplot(2,2,3); imagesc(img1); axis image;box off; hold on;
% plot(pupilStart1(1), pupilStart1(2),'r+')
% subplot(2,2,4); imagesc(img2); axis image; box off; hold on
% plot(pupilStart2(1), pupilStart2(2),'r+')

end

