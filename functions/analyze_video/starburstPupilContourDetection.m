% Starburst Algorithm
%
% This source code is part of the starburst algorithm.
% Starburst algorithm is free; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% Starburst algorithm is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with cvEyeTracker; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%
% Starburst Algorithm - Version 1.0.0
% Part of the openEyes ToolKit -- http://hcvl.hci.iastate.edu/openEyes
% Release Date:
% Authors : Dongheng Li <donghengli@gmail.com>
%           Derrick Parkhurst <derrick.parkhurst@hcvl.hci.iastate.edu>
% Copyright (c) 2005
% All Rights Reserved.

function [epx, epy, edge_thresh] = starburst_pupil_contour_detection(I, cx, cy, edge_thresh,radiiPupil,minfeatures)

% Input
% I = input image
% cx, cy = central start point of the feature detection process
% edge_thresh = best guess for the pupil contour threshold (30)

% Ouput
% epx = x coordinate of feature candidates [row vector]
% epy = y coordinate of feature candidates [row vector]

% edge_thresh = 10;               % edge_threshold = best guess for the pupil contour threshold (30)
min_edge_thresh = 2;
N = 50;                         % number of rays to use to detect feature points % 50
minimum_candidate_features=N*minfeatures;  % minimum number of pupil feature candidates ***6/10
dis = radiiPupil(1)*.65; % Distance from pupil center guess to start searching for edge
angle_spread = 60*pi/180; % Changed from 100*pi/180
loop_count = 0;
max_loop_count = 10;
tcx(loop_count+1) = cx;
tcy(loop_count+1) = cy;
angle_step= 2*pi/N;

while edge_thresh > min_edge_thresh && loop_count <= max_loop_count
    epx = [];
    epy = [];
    while length(epx) < minimum_candidate_features && edge_thresh > min_edge_thresh
        [epx, epy, epd] = locate_edge_points(I, cx, cy, dis, angle_step, 0, 2*pi, edge_thresh);
        if length(epx) < minimum_candidate_features
            edge_thresh = edge_thresh - 1;
        end
    end
    %     if edge_thresh <= min_edge_thresh
    %         break;
    %     end
    
    angle_normal = atan2(cy-epy, cx-epx);
    
    for i=1:length(epx)
        [tepx, tepy, ~] = locate_edge_points(I, epx(i), epy(i), dis*1.5, angle_step*(edge_thresh/epd(i)), angle_normal(i), angle_spread, edge_thresh);
        epx = [epx tepx];
        epy = [epy tepy];
    end
    
    loop_count = loop_count+1;
    tcx(loop_count+1) = mean(epx);
    tcy(loop_count+1) = mean(epy);
    
    if abs(tcx(loop_count+1)-cx) + abs(tcy(loop_count+1)-cy) < 10 % threshold?
        break;
    end
    cx = mean(epx);
    cy = mean(epy);
end

%% Plot final points used - debugging
%     figure(2); clf
%     imagesc(I);colormap(gray)
%     hold on;
%     plot([cx*ones(1,length(epx)); epx(:)'], [cy*ones(1,length(epy));epy(:)'],'m')
%     hold on; plot(cx,cy,'+w','LineWidth',3)
%     plot(epx,epy,'.y')
%     axis image

%% Warning messages
if loop_count > 10 %***
    fprintf('Warning! edge points did not converge in %d iterations.',loop_count);
    return;
end;

if edge_thresh <= min_edge_thresh %***
    fprintf('Warning! Adaptive threshold is too low!\n');
    return;
end


function [epx, epy, dir] = locate_edge_points(I, cx, cy, dis, angle_step, angle_normal, angle_spread, edge_thresh)

[height, width] = size(I);
epx = [];
epy = [];
dir = [];
ep_num = 0;  % ep stands for edge point
step = 3    ; %***2
maxDist = 2;
for angle=(angle_normal-angle_spread/2+0.0001):angle_step:(angle_normal+angle_spread/2)
   cosangle = cos(angle);
    sinangle = sin(angle);
    
    maxstep = round(dis*maxDist/step);
    dists = dis + step*(0:maxstep)';
    ps = [round(cx + dists*cosangle)  round(cy + dists*sinangle)]; %***
    
    mask =  ps(:,2) < height & ps(:,2) > 1 & ps(:,1) < width & ps(:,1) > 1;
    ps = ps(mask,:);
    
    vals = I(sub2ind([height, width], ps(:,2), ps(:,1)));    
    d = diff(vals);    
    ind = find(d>edge_thresh,1);
    
    if ~isempty(ind)
        ep_num = ep_num+1;
        epx(ep_num) = ps(ind,1)+step/2;    % edge point x coordinate
        epy(ep_num) = ps(ind,2)+step/2;    % edge point y coordinate
        dir(ep_num) = d(ind);        
    end    
    
end

%% DEBUGGING/VISUALIZATION %%%
%     figure(3); clf;
%     imagesc(I);colormap(gray); axis image
%     hold on; plot(epx, epy,'.y');
%     if angle_spread  <6;  c = 'c'; else  c = 'm'; end
% 
%     plot([cx*ones(1,length(epx)); epx], [cy*ones(1,length(epy));epy],c)
% 
%     hold on; plot(cx,cy,'om')
%     pause(.001)
