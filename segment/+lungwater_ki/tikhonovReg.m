function tikhonovReg(no, lambda)

global SET

if nargin ==1
    lambda=[];
end


lungMask=squeeze(SET(no).LungWater.LungMask);

for tf=1:SET(no).TSize
    image = squeeze(SET(no).IM(:,:,tf,:)); %scaled 0-1
   
    if isempty(lambda)
        midSlice=round(size(image,3)/2);
        lambda= lungwater_ki.LCurveFind(image(:,:,midSlice));
    end

    bodymask=zeros(size(image));
    bodymask(image>=mean(image(:)))=1;
    bodymask(squeeze(lungMask(:,:,tf,:))==1)=0; %exclude lungwater
    
    if tf==1 %same shading map for all timeframes to save processing times
        shading_map=zeros(size(image));
        for jj=1:size(image,3) %slice-by-slice Tikhonov regularization
            shading_map(:,:,jj) = lungwater_ki.tikReg2D(image(:,:,jj).*bodymask(:,:,jj),lambda);
        end
    end
    
    if isfield(SET(no).LungWater, 'ImageMOCO')
        tikhonovIm=squeeze(SET(no).LungWater.ImageMOCO(:,:,tf,:))./shading_map;
    else %normalize uncorrected image (DONT DO THIS FOR TIME-RESOLVED IMAGES)
        tikhonovIm=image./shading_map;
    end
    tikhonovIm(isinf(tikhonovIm))=0; tikhonovIm(isnan(tikhonovIm))=0;
    SET(no).LungWater.Lambda(1,tf)=lambda;
    SET(no).LungWater.ImageTikhonov(:,:,tf,:)=tikhonovIm;
end