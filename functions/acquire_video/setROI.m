function vid = setROI(vid, img1, img2,res, box)

%% Display snapshot
figure(1); 
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
h(1) = subplot(1,2,1);
imshow(img1); title('Camera 1')

h(2) = subplot(1,2,2);
imshow(img2); title('Camera 2')

%% Prompt for ROI around each eye
disp('Click on the camera coming from the left angle')
a = ginput(1);
hout = gca;

if hout == h(2) % right subplot clicked, switch cameras
    vid = [vid(2) vid(1)];
    img2temp = img1;
    img1 = img2;
    img2 = img2temp;
    subplot(1,2,1); imshow(img1); title('Camera 1')
    subplot(1,2,2); imshow(img2); title('Camera 2')    
end


disp('Drag a rectangle around the eye in camera 1')
subplot(1,2,1)
hRect = imrect;
pos1 = round(getPosition(hRect));

posnew(1) = pos1(2);
posnew(2) = pos1(1);
posnew(3) = pos1(4);
posnew(4) = pos1(3);


%% Set new ROI
roi = posnew;

roi(4) = 2*max(box/2-roi(2), roi(2)+roi(4)-box/2); % Ensure new image is still centered
roi(2) = box/2-roi(4)/2;

roi2 = roi + [(res-box)/2 0 0] + [200 0 0 0]; % vertical offset

set(vid(1), 'ROIPosition',roi2);
set(vid(2), 'ROIPosition',roi2);
close gcf


%% Write large images to file
imwrite(img1, 'img1large.png','png')
imwrite(img2, 'img2large.png','png')

