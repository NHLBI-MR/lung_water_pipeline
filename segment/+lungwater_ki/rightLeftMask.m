function [rightLungMask, leftLungMask]=rightLeftMask(no,tf, doPlot)
%divides lung mask into right and left lung

global SET

if nargin <2
    tf =1;
end

if nargin <3
    doPlot =0;
end


lungMask=squeeze(SET(no).LungWater.LungMask(:,:,tf,:));

cc = bwconncomp(lungMask,26);
labeled = labelmatrix(cc);

numPixels = cellfun(@numel,cc.PixelIdxList);
[biggest,idx] =sort(numPixels,'descend');

if cc.NumObjects==1
    slices=find(sum(sum(lungMask,1),2)~=0); %slices with mask
    if slices>1
        right_slices=slices(1:floor(length(slices)/2));
        left_slices=slices(find(slices==right_slices(end))+1:end);
    else
        right_slices=1;
        left_slices=1;
    end
    rightLungMask=zeros(size(lungMask));
    rightLungMask(:,:,right_slices)=lungMask(:,:,right_slices);

    leftLungMask=zeros(size(lungMask));
    leftLungMask(:,:,left_slices)=lungMask(:,:,left_slices);
    return;
end

S = regionprops(cc,'Centroid');
if size(S,1)>0
    pos1=S(idx(1)).Centroid(1);
    pos2=S(idx(2)).Centroid(1);

    if pos1(1)<pos2(1) % fist object is right lung
        rightLungMask=lungMask.*(labeled==idx(1));
        leftLungMask=lungMask.*(labeled==idx(2));
    else % fist object is left lung
        rightLungMask=lungMask.*(labeled==idx(2));
        leftLungMask=lungMask.*(labeled==idx(1));
    end
    if doPlot

        figure(89); clf;
        subplot(1,3,1); montage(labeled,'DisplayRange', [0, max([idx(1:2)])]);
        subplot(1,3,2); montage(rightLungMask,'DisplayRange', [0, 1])
        subplot(1,3,3); montage(leftLungMask,'DisplayRange', [0, 1])
    end
else
    rightLungMask=lungMask.*0;
    leftLungMask=lungMask.*0;
end