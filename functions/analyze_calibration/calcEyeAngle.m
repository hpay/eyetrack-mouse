function varargout = calcEyeAngle(results, theta)

%% Angle between two cameras in degrees
if ~exist('theta','var')
    theta = 40;
end

%% Horizontal

% Option 1: Original
d1 =  results.pupil1(:,1) - results.cr1a(:,1);
d2 = -(results.pupil2(:,1) - results.cr2b(:,1));

% Option 2: average 2 CRs
% d1 =  results.pupil1(:,1) - (results.cr1a(:,1)+results.cr1b(:,1))/2;
% d2 = -(results.pupil2(:,1) - (results.cr2a(:,1)+results.cr2b(:,1))/2);

% Normalize differences in magnification in two cameras?
% d1 = d1/nanmean(results.cr1b(:,1)-results.cr1a(:,1));
% d2 = d2/nanmean(results.cr2b(:,1)-results.cr2a(:,1));        

% Calculate angular position
posH = atand(d1./d2*sind(theta) ./ (1 + d1./d2*cosd(theta)));
% posH = atand(sind(theta) ./ (d2./d1 + cosd(theta))); % should be equivalent

% Desaccade mask for finding Rp more accurately
thresh = 5;
mask1 = filter(ones(10,1),1, abs([0; diff(d1)]) > thresh) > 0;
mask2= filter(ones(10,1),1, abs([0; diff(d2)]) > thresh) > 0;

rP = calcEyeRadius(d1, d2,theta, mask1|mask2)

% ALTERNATE METHOD - only needs one eye after rP calculation
posH1_alternate = (asind( d1/rP));
posH1_alternate(imag(posH1_alternate)~=0) = NaN;

posH2_alternate = (asind( -d2/rP)) + theta;
posH2_alternate(imag(posH2_alternate)~=0) = NaN;

% DEBUG - compare results using rP
figure; plot(posH); hold on; plot(posH1_alternate);
plot(posH2_alternate); 
legend('posH','posH1 using r_P','posH2 using r_P')

%% Vertical (less exact but pretty close)
if nargout>1
        width = 30; % width of medfilt for removing alignment points
    removeoutliers = @(x) (abs(x-medfilt1(x,width)) > 2*nanstd(medfilt1(x,width)));
    
    
    cr1y = results.cr1a(:,2);
    cr1y(removeoutliers(cr1y))= NaN;
    cr1y = inpaint_nans(cr1y);
    cr1y = medfilt1(cr1y,3);
    
    cr2y = results.cr2b(:,2);
    cr2y(removeoutliers(cr2y))= NaN;
    cr2y = inpaint_nans(cr2y);
    cr2y = medfilt1(cr2y,3);
    
    posV1 = asind((results.pupil1(:,2) - cr1y)./rP); % in degrees cam 1 estimate
    posV2 = asind((results.pupil2(:,2) - cr2y)./rP); % in degrees cam 2 estimate
    
    posV1(imag(posV1)~=0) = NaN;
    posV2(imag(posV2)~=0) = NaN;
    
    if nanstd(posV1) < nanstd(posV2)
        posV = posV1;
    else
        posV = posV2;
    end
end

%% Outputs
if nargout == 1
    varargout = {posH};
elseif nargout == 2
    varargout = {posH, posV};
else
    varargout = {posH, posV1, posV2};
end



end

function rP = calcEyeRadius(d1, d2, theta, mask) % Calculate pupil radius using approxiamate method of Stahl 2000
% rP_all = abs(d1-d2)/(2*sind(theta/2)); 
rP_all = mean([d1 d2],2)/sind(theta/2); 

% Remove desaccaded sections
if exist('mask','var')
    rP_all(mask)  = NaN;
end

% Remove outliers - 3 sigma
outliers = abs(rP_all-nanmedian(rP_all)) > 3*nanstd(rP_all);
rP_all(outliers) = NaN;

% Report final value
rP = nanmedian(rP_all);
end