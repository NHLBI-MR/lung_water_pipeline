function LWDMap_overlay(lwd_map,mask,image,titlestr,fignbr, noReshape, maxLim, minLim, tf)

if nargin<5
    fignbr=[];
end

if nargin<6
    noReshape=0;
end

if nargin<7
    maxLim=100;
end

if nargin<8
    minLim=0;
end


% slices=sort(unique(find(squeeze(nansum(nansum(squeeze(mask(:,:,tf,:)))))>0)))';  %slices with mask
if ndims(mask)>2
slices=1:size(mask,ndims(mask));
else
slices=1;
end

i = find(mask ==0);
mask(i) = NaN;


if ndims(mask)==4
    lwdMap=squeeze(lwd_map(:,:,tf,slices));
    im = squeeze(image(:,:,tf,slices));
    mask=squeeze(mask(:,:,tf,slices));
else
    lwdMap=squeeze(lwd_map(:,:,slices));
    im = squeeze(image(:,:,slices));
    mask=squeeze(mask(:,:,slices));
end
lwdMap= flipdim(lwdMap,3);
im= flipdim(im,3);
mask= flipdim(mask,3);

num=size(lwdMap,3);

%Transform RA images to one matrix
A  = im; %permute(repmat(1:size(im,3), [size(im,1) 1 size(im,3)]), [1 3 2]); %demo matrix
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
images = cell2mat(B); %convert to matrix

%Transform O_2 map into one mask matrix
A  = lwdMap.*mask; %permute(repmat(1:size(im,3), [size(im,1) 1 size(im,3)]), [1 3 2]); %demo matrix
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
lungWater = cell2mat(B); %convert to matrix

if~isempty(fignbr)
    figure(fignbr); clf;
else
    figure;
end

haxes=axes;
clim=[minLim,maxLim];
[hf,hb] = lungwater_ki.imoverlay2(images,lungWater,clim,[0,prctile(images(:),99)],'parula',1,haxes); %
if nargin>3
    title(titlestr,'FontSize',14, 'Color','k');
end
end