% Radial_Sym_Transform - Loy and Zelinski's fast radial feature detector
%
% An implementation of Loy and Zelinski's fast radial feature detector
%
% Usage: S = fastradial(im, radii, alpha);
%
% Arguments:
%            im    - image to be analysed (after using imread for example:im=imread('lenna.jpg');)
%            radii - array of integer radius values to be processed
%                    suggested radii might be [1 3 5]
%            alpha - radial strictness parameter.
%                    1 - slack, accepts features with bilateral symmetry.
%                    2 - a reasonable compromise.
%                    3 - strict, only accepts radial symmetry.
%                        ... and you can go higher
%
% Returns    S      - Symmetry map.  Bright points with high symmetry are
%                     marked with large positive values. Dark points of
%                     high symmetry marked with large negative values.
%
% To localize points use NONMAXSUPPTS on S, -S or abs(S) depending on
% what you are seeking to find.

% Reference:
% Loy, G.  Zelinsky, A.  Fast radial symmetry for detecting points of
% interest.  IEEE PAMI, Vol. 25, No. 8, August 2003. pp 959-973.

function S = radialSymTransform(im, radii, alpha, ploton)

if ~exist('ploton', 'var')
    ploton = 0;
end

    if any(radii ~= round(radii)) || any(radii < 1)
        error('radii must be integers and >= 1')
    end

    % Prevent dark spots detected near edges (HP 6/24/13)
%     maxval = NaN;%max(im(:));
%     im(1,:) = maxval; im(end,:) = maxval; im(:,1) = maxval; im(:,end) = maxval;

    
    [rows,cols]=size(im);

    % Use the Sobel masks to get gradients in x and y
    gx = [-1 0 1
          -2 0 2
          -1 0 1];
    gy = gx';

    imgx = filter2(gx,im);
    imgy = filter2(gy,im);
    mag = sqrt(imgx.^2 + imgy.^2)+eps; % (+eps to avoid division by 0)

    % Normalise gradient values so that [imgx imgy] form unit
    % direction vectors.
    imgx = imgx./mag;
    imgy = imgy./mag;

    S = zeros(rows,cols);  % Symmetry matrix

    [x,y] = meshgrid(1:cols, 1:rows);

    for n = radii
	M = zeros(rows,cols);  % Magnitude projection image
	O = zeros(rows,cols);  % Orientation projection image

        % Coordinates of 'positively' and 'negatively' affected pixels
        posx = x + round(n*imgx);
        posy = y + round(n*imgy);

        negx = x - round(n*imgx);
        negy = y - round(n*imgy);

        % Clamp coordinate values to range [1 rows 1 cols]
        posx( (posx<1) )    = 1;
        posx( (posx>cols) ) = cols;
        posy( (posy<1) )    = 1;
        posy( (posy>rows) ) = rows;

        negx( (negx<1) )    = 1;
        negx( (negx>cols) ) = cols;
        negy( (negy<1) )    = 1;
        negy( (negy>rows) ) = rows;


        % Form the orientation and magnitude projection matrices
        for r = 1:rows
            for c = 1:cols
                O(posy(r,c),posx(r,c)) = O(posy(r,c),posx(r,c)) + 1;
                O(negy(r,c),negx(r,c)) = O(negy(r,c),negx(r,c)) - 1;

                M(posy(r,c),posx(r,c)) = M(posy(r,c),posx(r,c)) + mag(r,c);
                M(negy(r,c),negx(r,c)) = M(negy(r,c),negx(r,c)) - mag(r,c);
            end
        end

        % Clamp Orientation projection matrix values to a maximum of
        % +/-kappa,  but first set the normalization parameter kappa to the
        % values suggested by Loy and Zelinski
        if n == 1, kappa = 8; else kappa = 9.9; end

        O(O >  kappa) =  kappa;
        O(O < -kappa) = -kappa;

        % Unsmoothed symmetry measure at this radius value
        F = M./kappa .* (abs(O)/kappa).^alpha;

        % Generate a Gaussian of size proportional to n to smooth and spread
        % the symmetry measure.  The Gaussian is also scaled in magnitude
        % by n so that large scales do not lose their relative weighting.
        A = fspecial('gaussian',[n n], 0.25*n) * n;

        S = S + filter2(A,F);

    end  % for each radius

    S = S/length(radii);  % Average

    if ploton
     figure
     subplot(2,1,1)
    colormap(gray);
    imagesc(im,[0 255]);
    title('input image');

     subplot(2,1,2)
    colormap(gray);
    imagesc(S);
    title(' Symmetry map');
    end
