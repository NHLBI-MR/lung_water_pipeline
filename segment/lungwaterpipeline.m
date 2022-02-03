function lungwaterpipeline(no)
%Plugin for Segment to run lung water pipeline developed by Felicia Seemann, 2022
%felicia.seemann@nih.gov

%no is the set number for the image intended for lung water analysis

global SET DATA NO

localPath='C:\Users\username\lung_water_pipeline\python'; %enter the local path to your neural network location

lung_segmentation(no, localPath);
liver_roi(no);
coil_shading_correction(no);
calc_lwd(no);

%display results
d=sprintf('Lung water density: %0.1f %%', SET(no).Developer.LWD);
LWD_overlay_maps(squeeze(SET(no).Developer.LWD_map),squeeze(SET(no).Developer.LungMaskML),squeeze(SET(no).IM),d,10, 0, 40);
end


function mask = lung_segmentation(no, localPath)
%this function calls a lung segmentation neural network in python through
%an external system command

global SET

%define paths


script=[localPath, 'compute_lung_seg.py'];
input=[localPath, 'image_ML.mat'];
output=[localPath, 'mask_ML.mat'];
command = sprintf('python3 %s %s %s', script,input,output);


image = squeeze(calcfunctions('calctruedata',SET(no).IM, no));
%downsample
im=imresize(image, round(SET(no).ResolutionX/1.5), 'bicubic');

%save image as .mat
disp('Saving image to server');
save([localPath, 'image_ML.mat'], 'im');

%run model in python
disp('Running neural network');
system(command);

%pause to make sure mask has time to be saved
pause(5);

%load generated mask
disp('Loading mask from server');
load([localPath, 'mask_ML.mat']);

%downsample
mask=imresize(mask, 1/round(SET(no).ResolutionX/1.5), 'bicubic');
mask(mask>=0.3)=1;
mask(mask<0.3)=0;

%save mask
SET(no).Developer.LungMaskML=[];
SET(no).Developer.LungMaskML(:,:,1,:)=mask;

%store in segment
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
image = squeeze(SET(no).IM); %scaled 0-1
lungMask=squeeze(SET(no).Developer.LungMaskML);
liverMask=squeeze(SET(no).Developer.LiverMask);
liverMaskSlice=find(sum(sum(liverMask,1),2)~=0);
 
%find mask of body, removing the lungs
bodymask=zeros(size(image));
bodymask(image>=mean(image(:)))=1;
bodymask(lungMask==1)=0; %exclude lungs

lambda = findLambdaSmoothing(image(:,:,liverMaskSlice));

coilmap=[];
for loop=1:size(image,3) %slice-by-slice Tikhonov regularization
    coilmap(:,:,loop) = tikReg2D(image(:,:,loop).*bodymask(:,:,loop),lambda);
end

%normalize image
 image_norm=image./coilmap;

 image_norm(find(isinf(image_norm)))=0;
 image_norm(find(isnan(image_norm)))=0;

SET(no).Developer.ImageNorm(:,:,1,:)=image_norm;
SET(no).Developer.Coilmap(:,:,1,:)=coilmap;
SET(no).Developer.Lambda=lambda;



figure(4); clf;
ax(1)=subplot(2,3,1);  montage(image, 'DisplayRange', [0,max(image(:))]); title('No normalization')
ax(2)=subplot(2,3,2);  montage(image_norm, 'DisplayRange', [0,max(image(:))]); title('With normalization')
ax(3)=subplot(2,3,3);  montage(coilmap, 'DisplayRange', [0,max(coilmap(:))]);  title('Coil map')

ax(4)=subplot(2,3,4);  montage(permute(image,[1 3 2]), 'DisplayRange', [0,max(image(:))]); title('No normalization')
ax(5)=subplot(2,3,5);  montage(permute(image_norm,[1 3 2]), 'DisplayRange', [0,max(image(:))]); title('With normalization')
ax(6)=subplot(2,3,6);  montage(permute(coilmap,[1 3 2]), 'DisplayRange', [0,max(coilmap(:))]);  title('Coil map')

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
catch %if normalized image no found, pull original image
    SET(no).Developer.ImageNorm = calcfunctions('calctruedata',SET(no).IM, no); 
end

im=squeeze(SET(no).Developer.ImageNorm);
makeLungMaskFromRoi(no,'Lung');
lungMask=squeeze(SET(no).Developer.LungMaskML);
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

function lungmask2roi(no)
global SET DATA NO

NO=no;
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

[mi,mj] = ind2sub(size(mask),find(mask,1));

X = bwtraceboundary(mask,[mi,mj],'W');

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



function lambda = findLambdaSmoothing(image, lambda_range)
%this function is an adaptation of the function LCurveFind impelemned by
% W. Quinn Meadus, June 2019, which is avaliable at 
% https://github.com/meadus/TikhonovRegularizationSurfaceFit

% LCurveFind() finds a smoothing parameter to balance a tikhonov regularization problem
% based on the L-curve method as described in the 1992 paper. This function
% was used to determine a proper smoothing value for the lung water image
% set normalization. Effectively finds a balacance between error from
% fitting to the original data and error from oversmoothing.
%
% [spf,ind,xL,yL] = LCurveFind(slice)
% slice == slice to be normalized, zero values are excluded from the fit
% spf == final smoothing parameter (lambda)
% ind === index of array of possible smoothing parameters
% xL, yL == x and y for the L-curve
% plot(xL,yL,xL(ind),yL(ind),'*') to see what point of the curve was chosen

if nargin < 2
   lambda_range=linspace(0.5,2000); %lambda search range
else
    lambda_range=linspace(lambda_range(1),lambda_range(2)); %lambda search range
end
image=double(image);

ind = image>0; %values less than zero will be excluded from the fitting process
b = image(image>0); 


for i = 1:length(lambda_range)
    %Fitting the surface with a specific smoothing parameter
    [x,A,T] = tikReg2D(image,lambda_range(i));

    %calculating the errors with that parameter
    res_norm(i) = norm(A*x(:)-b,'fro');
    solution_norm(i) = norm(T*x(:),'fro');
end

res_norm_log = log(res_norm);
solution_norm_log = log(solution_norm);

x_grid = 0.5:0.25:50000;
%interpolate norms
res_norm_log= spline(lambda_range,res_norm_log,x_grid);
solution_norm_log = spline(lambda_range,solution_norm_log,x_grid);

%calculating maximum curvature, derivatives
xL1 = gradient(res_norm_log);
yL1 = gradient(solution_norm_log);

xL2 = del2(res_norm_log);
yL2 = del2(solution_norm_log);

k = (xL2.*yL1-xL1.*yL2)./(xL1.^2+yL1.^2).^1.5; %curvature equations
[~,ind] = min(k);
lambda = x_grid(ind); %optimized lambda at max curvature

end

function [X,A,T] = tikReg2D(image,lambda)
% This function was implemented by W. Quinn Meadus, June 2019, and is
% avaliable at https://github.com/meadus/TikhonovRegularizationSurfaceFit
% 
% tikReg2D() generates a surface to fit to the data in "slice". Based on
% John D'Errico's Gridfit. Zeros are ignored from the fit, tikhonov
% regularization allows for fitting a surface over data with large holes or
% missing data.
%
% slice == data to fit a surface over (will fit to non-zero data in the
% array), 2D matrix
% lambda == smoothing paramter, the higher the value, the smoother the fit


image=double(image);
[ny, nx, nz] = size(image);

b = image(image(:)>0); %rhs data, assuming 0 values are to be excluded from the fit
bind = find(image(:)); %rhs location in full grid

nb = length(b);
ngrid = length(image(:));

%Holds the information for the location of each b value in the full grid
%(bInd) while having a row corresponding to each b value.
A = sparse((1:nb)',bind, ones(nb,1),nb,ngrid);


%difference approximation in y
[i,j] = meshgrid(1:nx,2:(ny-1));
ind = j(:) + ny*(i(:)-1);
len = length(ind);

T2 = sparse(repmat(ind,1,3), [ind-1,ind,ind+1], [-1*ones(len,1),2*ones(len,1),-1*ones(len,1)], ngrid,ngrid);


%difference approximation in x
[i,j] = meshgrid(2:(nx-1),1:ny);
ind = j(:) + ny*(i(:)-1);
len = length(ind);

T1 = sparse(repmat(ind,1,3), [ind-ny,ind,ind+ny], [-1*ones(len,1),2*ones(len,1),-1*ones(len,1)], ngrid,ngrid);

%Combining regularization (tikhonov) matrices
T = [T1;T2];

%appending zeros to the rhs
b = [b;zeros(size(T,1),1)];

%solving the minimization problem (tikhonov regularization solution)
AT = [A;lambda*T];
X = reshape((AT'*AT)\(AT'*b),ny,nx);
end
