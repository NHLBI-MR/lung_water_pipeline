function lungwaterpipeline(no)
global SET DATA NO


% lung_segmentation(no);
liver_roi(no);
% coil_shading_correction(no);
calc_lwd(no);

%display results
d=sprintf('Lung water density: %0.1f %%', SET(no).Developer.LWD);
lungwater.LWD_overlay_maps(squeeze(SET(no).Developer.LWD_map),squeeze(SET(no).Developer.LungMaskML),squeeze(SET(no).IM),d,10, 0, 40);


end


function mask = lung_segmentation(no)
%this function calls a lung segmentation neural network in python through
%an external system command

global SET

%define paths
windowsPath='\\hl-share.nhlbi.nih.gov\tmb\lab-campbell\_code_\python\dl_lung_segmentation\implementation\';
kauaiPath = '/mnt/lab-campbell/_code_/python/dl_lung_segmentation/implementation/';
script=[kauaiPath, 'compute_lung_seg.py'];
input=[kauaiPath, 'image_ML.mat'];
output=[kauaiPath, 'mask_ML.mat'];
command = sprintf('plink -ssh seemannfh@kauai -batch python3 %s %s %s', script,input,output);


image = squeeze(calcfunctions('calctruedata',SET(no).IM, no));
%downsample
im=imresize(image, round(SET(no).ResolutionX/1.25), 'bicubic');

%save image as .mat
disp('Saving image to server');
save([windowsPath, 'image_ML.mat'], 'im');

%run model in python
disp('Running neural network');
system(command);

%pause to make sure mask has time to be saved
pause(5);

%load generated mask
disp('Loading mask from server');
load([windowsPath, 'mask_ML.mat']);

%downsample
mask=imresize(mask, 1/round(SET(no).ResolutionX/1.25), 'bicubic');
mask(mask>=0.3)=1;
mask(mask<0.3)=0;

%save mask
SET(no).Developer.LungMaskML=[];
SET(no).Developer.LungMaskML(:,:,1,:)=mask;

%store in segment
% make_mask_overlay(image,mask);
lungmask2roi(no);
maskHoles2Roi(no);
makeLungMaskFromRoi(no,'Lung');

SET(no).Developer.ImageNorm =calcfunctions('calctruedata',SET(no).IM, no);
end

%---------------------------------------------------------------------
function liver_roi(no)
%---------------------------------------------------------------------
%this function places a circular ROI in the liver, under the right lung

global SET NO
mask=squeeze(SET(no).Developer.LungMaskML);
liverMask=zeros(size(mask));

% try
%0) shrink mask a bit for stability
dilatelevel=5;
for sl=1:size(mask,3)
    mask(:,:,sl) = imerode(mask(:,:,sl),strel('square',dilatelevel));
end

%1) find right lung
cc = bwconncomp(mask,26);
labeled = labelmatrix(cc);
numPixels = cellfun(@numel,cc.PixelIdxList);
[~,idx] =sort(numPixels,'descend');
S = regionprops(cc,'Centroid');
if size(S,1)>0
    pos1=S(idx(1)).Centroid(1);
    pos2=S(idx(2)).Centroid(1);
    
    if pos1(1)<pos2(1) % fist object is right lung
        rightLungMask=mask.*(labeled==idx(1));
        S = S(idx(1));
    else % fist object is left lung
        rightLungMask=mask.*(labeled==idx(2));
        S = S(idx(2));
    end
else
    rightLungMask=mask;
end

%2) find mid lung slice using centroid and middle column of right lung
mid=round(S(1).Centroid);
midSlice=mid(3);
midCol=mid(1);
bottomMidRow=find(sum(rightLungMask(:,:,midSlice)>0, 2)>0, 1, 'last');

%coordinates for center of ROI
r0=20; %radius for liver ROI, mm
n = 79; %number of points (segment uses 80 points)
omega = ((2*pi/n)*(1:n))';
rx = r0/SET(no).ResolutionX;
ry = r0/SET(no).ResolutionY;
X=bottomMidRow+ry+ceil(17.5/SET(no).ResolutionY); %place roi ~8.75mm below lung segmentation
Y=midCol;

%create circular mask
x = repmat(rx*sin(omega)+X,[1 1 1]); x = [x ; x(1,:)];
y = repmat(ry*cos(omega)+Y,[1 1 1]); y = [y ; y(1,:)];
liverMask(:,:,midSlice) = poly2mask(y,x,size(mask,1),size(mask,2));
% catch
%     liverMask(10:10,10:10,1) = 1;
% end

%delete any current liver ROI
roi('roitemplatedelete_Callback','Liver');

%store mask as ROI
NO=no;
SET(NO).CurrentTimeFrame=1;
n = size(SET(NO).Roi,2)+1;

SET(NO).Roi(n).X = x(:);
SET(NO).Roi(n).Y = y(:);
SET(NO).Roi(n).T = 1;
SET(NO).Roi(n).Z = midSlice;
SET(NO).Roi(n).Sign = 1;
[~,area] = calcfunctions('calcroiarea',NO,n);
SET(NO).Roi(n).Area = area;
[m,sd]=calcfunctions('calcroiintensity',NO,n);
SET(NO).Roi(n).Mean = m;
SET(NO).Roi(n).StD = sd;
SET(NO).Roi(n).Name = 'Liver';
SET(NO).Roi(n).LineSpec = 'r-';
SET(NO).RoiN = n;
SET(NO).RoiCurrent = 1;


SET(no).Developer.LiverMask=[];
SET(no).Developer.LiverMask(:,:,1,:)=liverMask;
end

%---------------------------------------------------------------------
function coil_shading_correction(no)
%---------------------------------------------------------------------
%this normalizes the images by fitting them to a 3rd order polynomial

global SET
image = squeeze(calcfunctions('calctruedata',SET(no).IM, no));
lungMask=squeeze(SET(no).Developer.LungMaskML);
liverMask=squeeze(SET(no).Developer.LiverMask);

%crop images around masks for stability
rows=find(sum(sum(lungMask+liverMask,3),2));
cols=find(sum(sum(lungMask+liverMask,3),1));
firstRow=rows(1)-3;  lastRow=rows(end)+10;
firstCol=cols(1)-10;  lastCol=cols(end)+10;

im = image(firstRow:lastRow,firstCol:lastCol,:);
lungMask=lungMask(firstRow:lastRow,firstCol:lastCol,:);

%find mask of body, removing the lungs
object_mask = im > mean(im(:)) - 1*std(im(:));
object_mask = object_mask - lungMask;
object_mask(object_mask<0) = 0;
object_mask(object_mask>0) = 1;

[coilmap, B]=polynomial_fit(im, object_mask); %3rd order polynomial fit over 3d volume

%normalize image
psn_im = im./coilmap;

image_norm=zeros(size(image));
image_norm(firstRow:lastRow,firstCol:lastCol,:) = psn_im;
SET(no).Developer.ImageNorm(:,:,1,:)=image_norm;
SET(no).Developer.Coilmap(:,:,1,:)=coilmap;
SET(no).Developer.NormCoeffs=B;


figure(4); clf;
ax(1)=subplot(2,3,1);  montage(image, 'DisplayRange', [0,max(im(:))]); title('No PSN')
ax(2)=subplot(2,3,2);  montage(image_norm, 'DisplayRange', [0,max(im(:))]); title('With PSN')
ax(3)=subplot(2,3,3);  montage(coilmap, 'DisplayRange', [0,max(coilmap(:))]);  title('Polynomial fit')

ax(4)=subplot(2,3,4);  montage(permute(image,[1 3 2]), 'DisplayRange', [0,max(im(:))]); title('No PSN')
ax(5)=subplot(2,3,5);  montage(permute(image_norm,[1 3 2]), 'DisplayRange', [0,max(im(:))]); title('With PSN')
ax(6)=subplot(2,3,6);  montage(permute(coilmap,[1 3 2]), 'DisplayRange', [0,max(coilmap(:))]);  title('Polynomial fit')

colormap(ax(1), 'gray');
colormap(ax(2),'gray');
colormap(ax(3),'parula');

colormap(ax(4), 'gray');
colormap(ax(5),'gray');
colormap(ax(6),'parula');

end

%---------------------------------------------------------------------
function calc_lwd(no)
%---------------------------------------------------------------------
%This function calculates the lung water density relative the liver

global SET
try
    SET(no).Developer.ImageNorm;
catch
    SET(no).Developer.ImageNorm =calcfunctions('calctruedata',SET(no).IM, no);
end



im=squeeze(SET(no).Developer.ImageNorm);
makeLungMaskFromRoi(no,'Lung');
% lungMask=squeeze(SET(no).Developer.LungMaskML);
lungMask=squeeze(SET(no).Developer.LungMask);
liverMask=squeeze(SET(no).Developer.LiverMask);

LWD_map = 70*im./mean(im(liverMask==1)); %assuming hepatic water density is 70%
LWD = mean(LWD_map(lungMask==1));
LWD_map((lungMask+liverMask)==0)=im((lungMask+liverMask)==0)/10;

SET(no).Developer.LWD_map=[]; SET(no).Developer.LWD=[];
SET(no).Developer.LWD_map(:,:,1,:)=LWD_map;
SET(no).Developer.LWD=LWD;
end

%---------------------------------------------------------------------
%SUPPORTING FUNCTIONS FOR SEGMENT
%---------------------------------------------------------------------

function [coilmap, B]=polynomial_fit(im, object_mask)

[Yp, Xp, Zp] = meshgrid(1:size(im,2), 1:size(im,1), 1:size(im,3));
masked_inds=find(object_mask);
x=Xp(masked_inds); y=Yp(masked_inds); z=Zp(masked_inds);
C=[x.^3, y.^3, z.^3, x.^2.*y, x.^2.*z, y.^2.*x, y.^2.*z, z.^2.*x, z.^2.*y, x.*y.*z, x.^2, y.^2, z.^2, x.*y, y.*z, x.*z, x, y, z, ones(size(x))];

si=double(im(masked_inds));

pn_factor       = 500; % normalised
search_range    = 10;
fit_minimum     = [repmat(-search_range, [1 19]) -1]*pn_factor;
fit_maximum     = [repmat(search_range, [1 19]) 1]*pn_factor;

opts=  optimset('Display','none');

B=lsqlin(C,si, [], [], [], [], fit_minimum, fit_maximum,[],opts)';

%polynomial in 3D coordinates
surfit_s = @(B,XYZ) B(1)*XYZ(:,:,:,1).^3 + B(2)*XYZ(:,:,:,2).^3 + B(3)*XYZ(:,:,:,3).^3 +                                        ... [3rd order] terms
    B(4)*XYZ(:,:,:,1).^2.*XYZ(:,:,:,2) + B(5)*XYZ(:,:,:,1).^2.*XYZ(:,:,:,3) +                                                   ... [3rd order] cross terms //x2
    B(6)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,1) + B(7)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,3) +                                                   ... [3rd order] cross terms //y2
    B(8)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,1) + B(9)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,2) + B(10)*XYZ(:,:,:,1).*XYZ(:,:,:,2).*XYZ(:,:,:,3) +  ... [3rd order] cross terms //z2
    B(11)*XYZ(:,:,:,1).^2 + B(12)*XYZ(:,:,:,2).^2 + B(13)*XYZ(:,:,:,3).^2 +                                                     ... [2nd order] terms
    B(14)*XYZ(:,:,:,1).*XYZ(:,:,:,2)  + B(15)*XYZ(:,:,:,2).*XYZ(:,:,:,3)  + B(16)*XYZ(:,:,:,1).*XYZ(:,:,:,3) +                  ... [2nd order] cross terms
    B(17)*XYZ(:,:,:,1) + B(18)*XYZ(:,:,:,2) + B(19)*XYZ(:,:,:,3) + B(20);

%3D coordinated for coil map
psn_phys_coords = zeros([size(im), 3], 'single');
psn_phys_coords(:,:,:,1) = Xp; psn_phys_coords(:,:,:,2) = Yp; psn_phys_coords(:,:,:,3) = Zp; %store meshgrid in 4D array

coilmap = surfit_s(B, psn_phys_coords);


cutoff_coilmap=prctile(coilmap(:),95); %cap at 95th percentile
norm_coilmap = coilmap./cutoff_coilmap;
lower_cutoff=0.1;
norm_coilmap(norm_coilmap<lower_cutoff) = lower_cutoff; %lowe cap at 0.1 or maybe higher
norm_coilmap(norm_coilmap>1)=1;
coilmap=norm_coilmap;

end


function lungmask2roi(no)
global SET DATA NO

NO=no;
% im=squeeze(SET(no).IM(:,:,1,:));

mask=squeeze(SET(no).Developer.LungMaskML);

DATA.Silent=true;
roi('roitemplatedelete_Callback','Lung');
DATA.Silent=false;

n=length(SET(no).Roi);
for slice=1:size(mask,3) %loop over slices
    cc = bwconncomp(mask(:,:,slice),8);
    labeled = labelmatrix(cc);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [biggest,idx] =sort(numPixels,'descend');
    inds=idx(biggest>5);
    
    for i=1:length(inds) %loop over binary objects in mask
        m=(labeled==inds(i));
        [x, y]=getObjectContours(m);
        
        if ~isempty(x)||~isempty(y)
            SET(NO).CurrentTimeFrame=1;
            %store objects and holes as ROIs
            n = n+1;
            SET(NO).Roi(n).X = x(:);
            SET(NO).Roi(n).Y = y(:);
            SET(NO).Roi(n).T = 1;
            SET(NO).Roi(n).Z = slice;
            SET(NO).Roi(n).Sign = 1;
            [~,area] = calcfunctions('calcroiarea',NO,n);
            SET(NO).Roi(n).Area = area;
            [m,sd]=calcfunctions('calcroiintensity',NO,n);
            SET(NO).Roi(n).Mean = m;
            SET(NO).Roi(n).StD = sd;
            SET(NO).Roi(n).Name = 'Lung';
            SET(NO).Roi(n).LineSpec = 'b-';
            SET(NO).RoiN = n;
            SET(NO).RoiCurrent = SET(NO).RoiN;
            %             roi('expandcontract_Callback',2);
        end
    end
end
SET(NO).RoiCurrent = SET(NO).RoiN;

panels = find(ismember(DATA.ViewPanels,SET(no).Linked));
for p = panels
    drawfunctions('drawroi',p);
end
end


function [x, y]=getObjectContours(mask)

global DATA
%finds the contours around the object (only) in mask
% mask = bwareafilt(logical(mask_in), 1);

[mi,mj] = ind2sub(size(mask),find(mask,1));

X = bwtraceboundary(mask,[mi,mj],'W');
%  x=[X(:,1),X(1,1)]'; y=[X(:,2),X(1,2)]'; %first point same as last

windowWidth = 45;
polynomialOrder = 12;

if size(X,1)<windowWidth
    if mod(size(X,1),2)~=0 %odd
        windowWidth = size(X,1)-2;
    else %even
        windowWidth = size(X,1)-1;
    end
end

if windowWidth>12
    x = sgolayfilt([X(end-ceil(windowWidth/2):end-1,1)',X(:,1)',X(2:ceil(windowWidth/2),1)'], polynomialOrder, windowWidth);
    y = sgolayfilt([X(end-ceil(windowWidth/2):end-1,2)',X(:,2)',X(2:ceil(windowWidth/2),2)'], polynomialOrder, windowWidth);
    x = x(ceil(windowWidth/2):end-ceil(windowWidth/2)+1);
    y = y(ceil(windowWidth/2):end-ceil(windowWidth/2)+1);
    [x,y] = calcfunctions('resamplecurve',x',y',DATA.NumPoints-1);
    x=[x,x(1)]'; y=[y,y(1)]'; %first point same as last
else
    x=[]; y=[];
end
end

%-----------------------------------------------------
function maskHoles2Roi(no)
%-----------------------------------------------------
%Finds holes in a mask and adds ROIs around them in SET(no)
global SET DATA NO

viewfunctions('switchimagestack',no);
tools('enableundo');

windowWidth=45;
polynomialOrder=12;

n = SET(NO).RoiN;
% figure(100);
for slice=1:size(SET(no).Developer.LungMaskML,4)
    BW = SET(no).Developer.LungMaskML(:,:,slice);
    [B,L,N] = bwboundaries(BW,'holes');
    
    for k=N+1:length(B)
        boundary = B{k};
        
        if size(boundary,1)>2
            try
                hull = convhull(boundary(:,1),boundary(:,2));
                x=boundary(hull,1);
                y=boundary(hull,2);
            catch
                x=boundary(:,1);
                y=boundary(:,2);
            end
            %             plot(y,x, 'w','LineWidth',2);
            x=[x; x(1)];
            y=[y; y(1)];
            [x,y] = calcfunctions('resampleclosedcurve',x,y,DATA.NumPoints-1);
            
            %store objects and holes as ROIs
            n = n+1;
            SET(NO).Roi(n).X = x';
            SET(NO).Roi(n).Y = y';
            SET(NO).Roi(n).T = 1;
            SET(NO).Roi(n).Z = slice;
            SET(NO).Roi(n).Sign = 1;
            [meanarea,area] = calcfunctions('calcroiarea',NO,n);
            SET(NO).Roi(n).Area = area;
            [m,sd]=calcfunctions('calcroiintensity',NO,n);
            SET(NO).Roi(n).Mean = m;
            SET(NO).Roi(n).StD = sd;
            SET(NO).Roi(n).Name = 'Lung exclude';
            SET(NO).Roi(n).LineSpec = 'w-';
            SET(NO).RoiN = n;
            SET(NO).RoiCurrent = SET(NO).RoiN;
            roi('expandcontract_Callback',2);
        end
    end
end
SET(NO).RoiN = n;
SET(NO).RoiCurrent = SET(NO).RoiN;
%Get all linked stacks and draw the rois in them as well.

panels = find(ismember(DATA.ViewPanels,SET(no).Linked));
for p = panels
    drawfunctions('drawroi',p);
end

end


function [roimask,lungVolume]=makeLungMaskFromRoi(no, str, mode)

global SET

if nargin ==1
    str='Lung';
end

if nargin < 3
    mode='Human';
end

mask=squeeze(SET(no).IM);
rois = indexROI(no,str); %lung ROI's
roimask=zeros(size(mask));

for rloop=rois
    if ~isnan(SET(no).Roi(rloop).Z)
        roimask(:,:,SET(no).Roi(rloop).Z) = roimask(:,:,SET(no).Roi(rloop).Z) + segment('createmask',...
            [size(mask,1) size(mask,2)],...
            SET(no).Roi(rloop).Y(:,1),...
            SET(no).Roi(rloop).X(:,1));
    end
end

roimask(roimask<=0)=0;
roimask(roimask>0)=1;

if strcmp(str,'Lung')
    SET(no).Developer.LungMask=[];
    SET(no).Developer.LungMask(:,:,1,:) = roimask;
elseif strcmp(str,'Vial100')
    SET(no).Developer.VialMask100=[];
    SET(no).Developer.VialMask100(:,:,1,:) = roimask;
elseif strcmp(str,'Vial50')
    SET(no).Developer.VialMask50=[];
    SET(no).Developer.VialMask50(:,:,1,:) = roimask;
elseif strcmp(str,'Liver')
    SET(no).Developer.LiverMask=[];
    SET(no).Developer.LiverMask(:,:,1,:) = roimask;
elseif strcmp(str,'LiverML')
    SET(no).Developer.LiverMaskML=[];
    SET(no).Developer.LiverMaskML(:,:,1,:) = roimask;
end


if strcmp(str,'Lung')
    lungVolume=sum(roimask(:))*SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap)*0.001; %in ml
else
    lungVolume=NaN;
end

end


%----------------------------------------------------------------
function [ind] = indexROI(no,str)
%----------------------------------------------------------------
%Find index of all ROI's named str
global SET
ind=[];
for i = 1:length(SET(no).Roi)
    if isequal(SET(no).Roi(i).Name,str)
        ind =[ind i];
    end
end

end


