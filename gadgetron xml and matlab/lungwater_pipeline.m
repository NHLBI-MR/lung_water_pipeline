function lungwater_pipeline(connection)
%this code performs lung water density analysis of a UTE image being
%streamed through gadgetron
%pathname is the path to yout local folder where the python network is stored
% data is a strcut containing image and image analysis

%Felicia Seemann 2022
%felicia.seemann@nih.gov

% !! USER ATTENTION !! Enter path to network!!
pathname = '.../lung_water_pipeline/python'; %local path of the python lung segmentation network

next_acquisition = @connection.next; %listens for images being streamed
data = pull_image(next_acquisition, connection); %pulles the binned image to analyze
data = lung_segmentation(data, pathname); %performs lung segmentation in python
data = liver_roi(data); %places an ROI in the liver
data = coil_shading_correction(data); %spatial normalization using tikhonov regularization
data = calc_lwd(data); %calculates the lung water density relative the liver

data=stream_results(data, connection); %streams back data to gadgetron

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUPPORTING FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
function [data] = pull_image(input, connection)
%pulles the binned image to analyze
%-------------------------------------------------------------------------
for i=1:2
    acquisition = input();
    acquisition.data=1000*(acquisition.data./max(acquisition.data(:))); %scale up low SIs coming from the scanner
    if i==2 %pull second image which is the binned image
        data.image = rot90(squeeze(acquisition.data),2);
        data.image_norm = data.image;
        data.reference = acquisition.header;
    end
    connection.send(acquisition);
end

end


%-------------------------------------------------------------------------
function data = lung_segmentation(data,pathname)
%performs lung segmentation in python - remember to set correct path
%-------------------------------------------------------------------------
script=[pathname, 'compute_lung_seg.py'];
input=[pathname, 'image_ML.mat'];
output=[pathname, 'mask_ML.mat'];
command = sprintf('python3 %s %s %s', script,input,output);

resolution=double(data.reference.field_of_view)./double(data.reference.matrix_size);

im=data.image;
im=imresize(im, round(resolution(1)/1.5), 'bicubic'); %upsample to ~1.5mm resolution
save(input, 'im');

system(command); %run network
load(output);

mask=imresize(mask, 1/round(resolution(1)/1.5), 'bicubic'); %downsample
mask(mask>=0.3)=1;
mask(mask<0.3)=0;

data.lungMask=mask;

end

%-------------------------------------------------------------------------
function data = liver_roi(data)
%places an ROI in the liver
%-------------------------------------------------------------------------
mask=data.lungMask;
liverMask=zeros(size(mask));

try
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
    resolution=double(data.reference.field_of_view)./double(data.reference.matrix_size);
    r0=20; %radius for liver ROI, mm
    n = 79; %number of points (segment uses 80 points)
    omega = ((2*pi/n)*(1:n))';
    rx = r0/resolution(1);
    ry = r0/resolution(2);
    X=bottomMidRow+ry+ceil(17.5/resolution(2)); %place roi ~8.75mm below lung segmentation
    Y=midCol;

    %create circular mask
    x = repmat(rx*sin(omega)+X,[1 1 1]); x = [x ; x(1,:)];
    y = repmat(ry*cos(omega)+Y,[1 1 1]); y = [y ; y(1,:)];
    liverMask(:,:,midSlice) = poly2mask(y,x,size(mask,1),size(mask,2));
catch
    liverMask(10:10,10:10,1) = 1;
end
data.liverMask=single(liverMask);
end

%-------------------------------------------------------------------------
function data = coil_shading_correction(data)
%spatial normalization using tikhonov regularization
%-------------------------------------------------------------------------
liverMaskSlice=find(sum(sum(data.liverMask,1),2)~=0); %find slice with liver ROI

bodymask=zeros(size(data.image));
bodymask(data.image>=mean(data.image(:)))=1;
bodymask(data.lungMask==1)=0; %exclude lungs

im = (data.image-min(data.image(:)))./max((data.image-min(data.image(:)))); %normalize for 0-1


lambda= findLambdaSmoothing(im(:,:,liverMaskSlice));
% lambda=40.75;
shading_map=[];
for loop=1:size(im,3) %slice-by-slice Tikhonov regularization
    shading_map(:,:,loop) = tikReg2D(im(:,:,loop).*bodymask(:,:,loop),lambda);
end

%normalize image
im_normalized = im./shading_map;

 image_normalized(find(isinf(image_normalized)))=0;
 image_normalized(find(isnan(image_normalized)))=0;

data.image_norm=im_normalized;
data.shading_map=shading_map;
data.lambda=lambda;

end

%-------------------------------------------------------------------------
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
%-------------------------------------------------------------------------
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

%-------------------------------------------------------------------------
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
%-------------------------------------------------------------------------
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

%-------------------------------------------------------------------------
function data = calc_lwd(data)
%calculates the lung water density relative the liver
%-------------------------------------------------------------------------
im=data.image_norm;

data.LWD_map = 70*im./mean(im(data.liverMask==1)); %assuming hepatic water density is 70%
data.LWD = mean(data.LWD_map(data.lungMask==1));
data.LWD_map((data.lungMask+data.liverMask)==0)=im((data.lungMask+data.liverMask)==0)/10;


end

%-------------------------------------------------------------------------
function data=stream_results(data, connection)
%streams back data to gadgetron
%-------------------------------------------------------------------------
%rotate back images and masks 
k=10;
data.LWD_map = k*rot90(single(round(data.LWD_map)),2);
data.masks=k*rot90(single(round(100*data.lungMask+70*data.liverMask)),2);
data.image = rot90(single(round(data.image)),2);
data.image_norm = rot90(single(round(data.image_norm)),2);

%stream back LWD map to gadgetron
header = connection.header;
image = gadgetron.types.Image.from_data(permute(data.LWD_map,[4 1 2 3]), data.reference);
image.header.field_of_view = [header.encoding.reconSpace.fieldOfView_mm.x header.encoding.reconSpace.fieldOfView_mm.y header.encoding.reconSpace.fieldOfView_mm.z];
image.header.image_type     = gadgetron.types.Image.MAGNITUDE;
image.header.repetition     = 0; 
image.header.image_index    = 1; 
image.header.image_series_index = 3;
connection.send(image);

%stream back lung and liver mask to gadgetron
header = connection.header;
image = gadgetron.types.Image.from_data(permute(data.masks,[4 1 2 3]), data.reference);
image.header.field_of_view = [header.encoding.reconSpace.fieldOfView_mm.x header.encoding.reconSpace.fieldOfView_mm.y header.encoding.reconSpace.fieldOfView_mm.z];
image.header.image_type     = gadgetron.types.Image.MAGNITUDE;
image.header.repetition     = 0; 
image.header.image_index    = 1; 
image.header.image_series_index = 4;
connection.send(image);

%stream back normalized image to gadgetron
header = connection.header;
image = gadgetron.types.Image.from_data(permute(data.image_norm,[4 1 2 3]), data.reference);
image.header.field_of_view = [header.encoding.reconSpace.fieldOfView_mm.x header.encoding.reconSpace.fieldOfView_mm.y header.encoding.reconSpace.fieldOfView_mm.z];
image.header.image_type     = gadgetron.types.Image.MAGNITUDE;
image.header.repetition     = 0; 
image.header.image_index    = 1; 
image.header.image_series_index = 5;
connection.send(image);
end
