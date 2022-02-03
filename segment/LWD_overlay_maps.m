function LWD_overlay_maps(im_overlay,mask,im_base,titlestr,fignbr, noReshape, maxLim);

%This function plots the image im_base with an overlay of im_overlay at the
%pixels indicated in mask. 

%Felicia Seemann 2022
%felicia.seemann@nih.gov

if nargin<5
    fignbr=100;
end

if nargin<6
    noReshape=0;
end

if nargin<7
    maxLim=100;
end
i = find(mask == 0); 
mask(i) = NaN; 

num=size(im_overlay,3); 


im_overlay= flipdim(im_overlay,3); 
mask= flipdim(mask,3); 
im_base= flipdim(im_base,3); 

%Transform base images to one matrix
A  = im_base;  %demo matrix
B = num2cell(A, [1 2]); %split the stack into a cell array 
if noReshape
    B = reshape(B, size(A,3)/num, num);
else
    if mod(size(A,3),4)==0
        B = reshape(B, 4, ceil(size(A,3)/4)); %reshape the tiles into their final position. You may need to transpose the reshape
    elseif mod(size(A,3),3)==0
        B = reshape(B, 3, ceil(size(A,3)/3));
    elseif mod(size(A,3),2)==0
        B = reshape(B, 2, ceil(size(A,3)/2));
    else
        B = reshape(B, size(A,3)/num, num);
    end
end
B1 = cell2mat(B); %convert to matrix

%Transform overlay image into one masked matrix
A  = im_overlay.*mask; %
B = num2cell(A, [1 2]); %split the stack into a cell array 
if noReshape
    B = reshape(B, size(A,3)/num, num); %reshape the tiles into their final position. You may need to transpose the reshape
else
    if mod(size(A,3),4)==0
        B = reshape(B, 4, ceil(size(A,3)/4)); %reshape the tiles into their final position. You may need to transpose the reshape
    elseif mod(size(A,3),3)==0
        B = reshape(B, 3, ceil(size(A,3)/3));
    elseif mod(size(A,3),2)==0
        B = reshape(B, 2, ceil(size(A,3)/2));
    else
        B = reshape(B, size(A,3)/num, num);
    end
end
B2 = cell2mat(B); %convert to matrix

figure(fignbr); clf; haxes=axes;
clim=[-5,maxLim];
[hf,hb] = imoverlay2(B1,B2,clim,[0,prctile(B1(:),99)],'parula',1,haxes); %
if nargin>3
    title(titlestr,'FontSize',14, 'Color','k');
end
end

function [hF,hB] = imoverlay2(B,F,climF,climB,cmap,alpha,haxes)
% IMOVERLAY(B,F) displays the image F transparently over the image B.
%    If the image sizes are unequal, image F will be scaled to the aspect
%    ratio of B.
%
%    [hF,hB] = imoverlay(B,F,[low,high]) limits the displayed range of data
%    values in F. These values map to the full range of values in the
%    current colormap.
%
%    [hF,hB] = imoverlay(B,F,[],[low,high]) limits the displayed range of
%    data values in B.
%
%    [hF,hB] = imoverlay(B,F,[],[],map) applies the colormap to the figure.
%    This can be an array of color values or a preset MATLAB colormaps
%    (e.g. 'jet' or 'hot').
%
%    [hF,hB] = imoverlay(B,F,[],[],[],alpha) sets the transparency level to
%    alpha with the range 0.0 <= alpha <= 1.0, where 0.0 is fully
%    transparent and 1.0 is fully opaque.
%
%    [hF,hB] = imoverlay(B,F,[],[],[],[],ha) displays the overlay in the
%    axes with handle ha.
%
%    [hF,hB] = imoverlay(...) returns the handles to the front and back
%    images.
%
%
% Author: Matthew Smith / University of Wisconsin / Department of Radiology
% Date created:  February 6, 2013
% Last modified: Jan 2, 2015
%
%
%  Examples:
%
%     % Overlay one image transparently onto another
%     imB = phantom(256);                       % Background image
%     imF = rgb2gray(imread('ngc6543a.jpg'));   % Foreground image
%     [hf,hb] = imoverlay(imB,imF,[40,180],[0,0.6],'jet',0.6);
%     colormap('parula'); % figure colormap still applies
%
%
%     % Use the interface for flexibility
%     imoverlay_tool;
%
%
% See also IMOVERLAY_TOOL, IMAGESC, HOLD, CAXIS.



ALPHADEFAULT = 0.4; % Default transparency value
CMAPDEFAULT = 'parula';

if nargin == 0,
    try
        imoverlay_tool;
        return;
    catch
        errordlg('Cannot find imoverlay_tool.', 'Error');
    end
end


% Check image sizes
if size(B,3) > 1
    error('Back image has %d dimensions!\n',length(size(B)));
end
if size(F,3) > 1
    error('Front image has %d dimensions!\n',length(size(F)));
end
if ~isequal(size(B),size(F))
    fprintf('Warning! Image sizes unequal. Undesired scaling may occur.\n');
end

% Check arguments
if nargin < 7
    haxes = [];
end

if nargin < 6 || isempty(alpha)
    alpha = ALPHADEFAULT;
end

if nargin < 5 || isempty(cmap)
    cmap = CMAPDEFAULT;
end

if nargin < 4 || isempty(climB)
    climB = [min(B(:)), max(B(:))];
end

if nargin < 3 || isempty(climF)
    climF = [min(F(:)), max(F(:))];
end

if abs(alpha) > 1
    error('Alpha must be between 0.0 and 1.0!');
end


% Create a figure unless axes is provided
if isempty(haxes) || ~ishandle(haxes)
    f=figure('Visible','off',...
        'Units','pixels','Renderer','opengl');
    pos = get(f,'Position');
    set(f,'Position',[pos(1),pos(2),size(B,2),size(B,1)]);
    haxes = axes;
    set(haxes,'Position',[0,0,1,1]);
    movegui(f,'center');
end
% Create colormap
cmapSize = 100;  % default size of 60 shows visible discretization
if ischar(cmap)
    
    try
        cmap = eval([cmap '(' num2str(cmapSize) ');']);
    catch
        fprintf('Colormap ''%s'' is not supported. Using ''jet''.\n',cmapName);
        cmap = jet(cmapSize);
    end
end
colormap(cmap);


% To have a grayscale background, replicate image to 3-channels
B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);

% Display the back image
axes(haxes);
%hB = imagesc(B);axis image off;
hB = imshow(B);axis image off;
% set(gca,'Position',[0,0,1,1]);

% Add the front image on top of the back image
hold on;
hF = imagesc(F,climF);

% If images are different sizes, map the front image to back coordinates
set(hF,'XData',get(hB,'XData'),...
    'YData',get(hB,'YData'))

% Make the foreground image transparent
alphadata = alpha.*(F >= climF(1));
set(hF,'AlphaData',alphadata);
if any(climF==0)
   colorbar('XTick',[climF(1), climF(2)],'FontSize',14, 'Color','k') 
else
    colorbar('XTick',[climF(1), 0, climF(2)],'FontSize',14, 'Color','k')
end
set(gcf,'color','w');

if exist('f')
    set(f,'Visible','on');
end


% Novel colormaps
%
% JET2 is the same as jet but with black base
    function J = jet2(m)
        if nargin < 1
            m = size(get(gcf,'colormap'),1);
        end
        J = jet(m); J(1,:) = [0,0,0];
    end

% JET3 is the same as jet but with white base
    function J = jet3(m)
        if nargin < 1
            m = size(get(gcf,'colormap'),1);
        end
        J = jet(m); J(1,:) = [1,1,1];
    end

% PARULA2 is the same as parula but with black base
    function J = parula2(m)
        if nargin < 1
            m = size(get(gcf,'colormap'),1);
        end
        J = parula(m); J(1,:) = [0,0,0];
    end

% HSV2 is the same as HSV but with black base
    function map = hsv2(m)
        map =hsv;
        map(1,:) = [0,0,0];
    end

% HSV3 is the same as HSV but with white base
    function map = hsv3(m)
        map =hsv;
        map(1,:) = [1,1,1];
    end

% HSV4 a slight modification of hsv (Hue-saturation-value color map)
    function map = hsv4(m)
        if nargin < 1, m = size(get(gcf,'colormap'),1); end
        h = (0:m-1)'/max(m,1);
        if isempty(h)
            map = [];
        else
            map = hsv2rgb([h h ones(m,1)]);
        end
    end
end

