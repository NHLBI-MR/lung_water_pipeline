function [roimask]=makeLungMaskFromRoiTimeResolved(no, str)

global SET
roimask=zeros(size(SET(no).IM));

if nargin<2
    str='Lung';
end

for slice=1:SET(no).ZSize
    rois = lungwater_ki.indexROISlice(no,str,slice); %find ROI's

    if ~isempty(rois) %loop over ROI's to make a mask
        rtemp=zeros(size(SET(no).IM,1), size(SET(no).IM,2), SET(no).TSize);
        for rloop=rois
            for tf=1:SET(no).TSize
                tf_tmp=tf;
                rtemp(:,:,tf)=rtemp(:,:,tf)+segment('createmask', [size(SET(no).IM,1) size(SET(no).IM,2)],SET(no).Roi(rloop).Y(:,tf_tmp),SET(no).Roi(rloop).X(:,tf_tmp));
            end
        end
        roimask(:,:,:, SET(no).Roi(rloop).Z) = roimask(:,:,:, SET(no).Roi(rloop).Z) + rtemp;
    end
end

%make sure mask is binary
roimask(roimask<=0)=0;
roimask(roimask>0)=1;


if isequal(str, 'Lung') %if lung mode, store mask to SET struct
    SET(no).LungWater.LungMask=[];
    SET(no).LungWater.LungMask = roimask;
    SET(no).LungWater.LungVolume=squeeze(sum(sum(sum(roimask(:,:,:,:))),4))'*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-6; %in Liter

    for tf=1:SET(no).TSize
        [SET(no).LungWater.RightLungMask(:,:,tf,:),SET(no).LungWater.LeftLungMask(:,:,tf,:)]=lungwater_ki.rightLeftMask(no,tf,0);
    end

end


