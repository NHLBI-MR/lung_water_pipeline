function LWD_overlay_maps(im_in,mask_in,im_RA,titlestr,fignbr, noReshape, maxLim);

if nargin<5
    fignbr=100;
end

if nargin<6
    noReshape=0;
end

if nargin<7
    maxLim=100;
end
i = find(mask_in ==0); 
mask_in(i) = NaN; 

num=size(im_in,3); 

im=im_in;
mask=mask_in;
imRA = im_RA;

im= flipdim(im,3); 
mask= flipdim(mask,3); 
imRA= flipdim(imRA,3); 

%Transform RA images to one matrix
A  = imRA; %permute(repmat(1:size(im,3), [size(im,1) 1 size(im,3)]), [1 3 2]); %demo matrix
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

%Transform O_2 map into one mask matrix
A  = im.*mask; %permute(repmat(1:size(im,3), [size(im,1) 1 size(im,3)]), [1 3 2]); %demo matrix
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
clim=[-5,50];
[hf,hb] = imoverlay2(B1,B2,clim,[0,prctile(B1(:),99)],'parula',1,haxes); %
if nargin>3
    title(titlestr,'FontSize',14, 'Color','k');
end
end
