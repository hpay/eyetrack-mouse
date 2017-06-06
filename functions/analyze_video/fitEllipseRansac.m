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

function [max_ellipse, max_inlier_indices, ransac_iter] = fitEllipseRansac(x, y, pupil_radii)

% Input
% [x,y] = feature points (row vectors)
% pupil_radius range of acceptable pupil
%
% Output
% max_ellipse = best fitting ellipse parameters
% max_inlier_indices = inlier indices for max_ellipse

max_ransac_iterations = 10000;     % maximum number of ransac iterations 10000


max_inliers = 0;
max_inlier_indices = [];
max_ellipse = [];
N = inf;
ransac_iter = 0;
nelements = 5;

x = x(:)';
y = y(:)';

if (isempty(x) || isempty(y))
    return;
end;

ep_num = length(x);
if (ep_num<5)
    fprintf(1,'Error! Must have at least 5 feature points to fit an ellipse\n');
    return
end;

[nx, ny, Hn] = normalize_point_coordinates(x, y);
dist_thres = 4*Hn(1); %Edit HP sqrt(3.84) dist_thres = Hn(1);
random_or_adaptive = 0;

ep = [nx; ny; ones(1,length(nx))]';

while N > ransac_iter || ransac_iter < 100 %% Added && ransac_iter < 100 HP 12/28/13
    if random_or_adaptive == 0
        % To ensure that 5 indices are different
        while 1
            random_indices = ceil(rand(nelements,1)*ep_num);
            rim = repmat(random_indices,[1 nelements]);
            if (sum(sum(rim==rim'))==nelements)
                break;
            end;
        end;
        % Calculate the ellipse from random selected points
        nxi=nx(random_indices);
        nyi=ny(random_indices);
    else
        % Calculate the ellipse from inliers of previous iteration
        nxi=nx(max_inlier_indices);
        nyi=ny(max_inlier_indices);
    end
    A = [nxi.*nxi; nxi.*nyi; nyi.*nyi; nxi; nyi; ones(1,length(nxi))]';
    
    [~, ~, va] = svd(A);
    nconic_par = va(:,end);
    nconic_matrix = [nconic_par(1) nconic_par(2)/2 nconic_par(4)/2;
        nconic_par(2)/2 nconic_par(3) nconic_par(5)/2;
        nconic_par(4)/2 nconic_par(5)/2 nconic_par(6)];
    diserr = sum((ep*nconic_matrix).*ep, 2);
    inliers_index = find(abs(diserr) < dist_thres);
    ninliers = length(inliers_index);
    
    random_or_adaptive = 0;
    if ninliers > max_inliers
        nellipse_par = convert_conic_parameters_to_ellipse_parameters(nconic_par);
        if isreal(nellipse_par) && nellipse_par(1) > 0 && nellipse_par(2) > 0
            ellipse_par = denormalize_ellipse_parameters(nellipse_par,Hn);
            er = min(ellipse_par(1:2))/max(ellipse_par(1:2));
            
            radius = max(ellipse_par(1:2));
            
            if  er > .8 % eccentricity of the ellipse % Changed from 0.75
                if  (radius>pupil_radii(1) && radius<pupil_radii(2)) %HP - check for correct sized pupil
                    
                    max_inliers = ninliers;
                    max_inlier_indices = inliers_index;
                    
                    % Edit to make a more intelligable structure: HP
                    temp = num2cell(ellipse_par);
                    [max_ellipse.a, max_ellipse.b, max_ellipse.x0, max_ellipse.y0, max_ellipse.beta] = deal(temp{:});
                    
                    N = log(1-0.99)/log(1-(ninliers/ep_num)^5+eps); % plus eps, in order to avoid log(0)
                    random_or_adaptive = 1;
                end
            end
        end
    end
    
    ransac_iter = ransac_iter+1;
    if (ransac_iter > max_ransac_iterations)
        error('Attention! The maximum number of ransac iterations exceeded!');
    end
end



end

function [nx, ny, H] = normalize_point_coordinates(x, y)

% Input:
% [x y] = coordinates (row vectors)

% Output:
% [nx ny] = normalized coordinates (row vectors)
% H = normalization homography

cx = mean(x);
cy = mean(y);
mean_dist = mean(sqrt(x.^2 + y.^2));
dist_scale = sqrt(2)/mean_dist;
H = [dist_scale     0            -dist_scale*cx;
    0              dist_scale   -dist_scale*cy;
    0              0            1];
nx = H(1,1)*x + H(1,3);
ny = H(2,2)*y + H(2,3);

end


function [e] = denormalize_ellipse_parameters(ne,Hn)

% Input
% ne = ellipse parameters constructed using normalized points
% Hn = homography used to normalize the points

% Output
% e = denormalized ellipse parameters

e(1) = ne(1) / Hn(1,1);
e(2) = ne(2) / Hn(2,2);
e(3) = (ne(3) - Hn(1,3)) / Hn(1,1);
e(4) = (ne(4) - Hn(2,3)) / Hn(2,2);
e(5)=ne(5);

end



function e = convert_conic_parameters_to_ellipse_parameters(c)

% This function converts conic parameters [c(1) c(2) c(3) c(4) c(5) c(6)
% to ellipse parameters [a b cx cy theta]

% Get ellipse orientation
theta = atan2(c(2),c(1)-c(3))/2;

% Get scaled major/minor axes
ct = cos(theta);
st = sin(theta);
ap = c(1)*ct*ct + c(2)*ct*st + c(3)*st*st;
cp = c(1)*st*st - c(2)*ct*st + c(3)*ct*ct;

% Get translations
T = [[c(1) c(2)/2]' [c(2)/2 c(3)]'];
t = -(2*T)\[c(4) c(5)]';
cx = t(1);
cy = t(2);

% Get scale factor
val = t'*T*t;
scale_inv = (val- c(6));

% Get major/minor axis radii
a = sqrt(scale_inv/ap);
b = sqrt(scale_inv/cp);

e = [a b cx cy theta];


if ~isreal(e)
    e = [0 0 0 0 0];
    return;
end
end