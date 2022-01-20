function lungwater_pipeline(connection)
disp("handle_connection was called.")


pathname = '/lung_water_pipeline/python'; %local path of the python lung segmentation code

next_acquisition = @connection.next;
data = pull_image(next_acquisition, connection);
data = lung_segmentation(data, pathname);
data = liver_roi(data);
data = calc_lwd(data);

% data
% figure(2); clf; montage(data.LWD_map, 'DisplayRange', [0,100]); colormap('parula');

data=stream_results(data, connection);

end

function [data] = pull_image(input, connection)

for i=1:2
    acquisition = input();
    acquisition.data=1000*(acquisition.data./max(acquisition.data(:))); %scale up low SIs coming from the scanner
    if i==2 %pull second image which is
        data.image = rot90(squeeze(acquisition.data),2);
        data.image_norm = data.image;
        data.reference = acquisition.header;
    end
    connection.send(acquisition);
end

end


function data = lung_segmentation(data,pathname)

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

function data = liver_roi(data)
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

function data = calc_lwd(data)
im=data.image_norm;

data.LWD_map = 70*im./mean(im(data.liverMask==1)); %assuming hepatic water density is 70%
data.LWD = mean(data.LWD_map(data.lungMask==1));
data.LWD_map((data.lungMask+data.liverMask)==0)=im((data.lungMask+data.liverMask)==0)/10;


end

function data=stream_results(data, connection)
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

end
