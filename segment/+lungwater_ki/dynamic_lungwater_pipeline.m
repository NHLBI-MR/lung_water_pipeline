function dynamic_lungwater_pipeline(nos, doLungSeg, doPlot, bodymaskslices, doLiverRoi, filename)
global SET NO DATA

disp('Post-processing pipeline');

if nargin<2
    doLungSeg=0;
end

if nargin<3
    doPlot=0;
end

if nargin<4
    bodymaskslices=[];
end

if nargin<5
    doLiverRoi=0;
end

if nargin<6
    filename='';
end



nbr_of_regions=10; %calculate anterior=posterior gradient based of 10 regions

for no=nos
    tf=1;
    try
        SET(no).LungWater.TimeVector=SET(no).LungWater.TimeVector(1:SET(no).TSize);
    catch
        SET(no).LungWater.TimeVector=SET(no).TimeVector;
    end
    disp('----------------------------------------');
    if doLungSeg %run neural network
        disp('Lung segmentation');
        SET(no).LungWater.LungMask=[];
        SET(no).LungWater.BodyMask=[];
        % [mask,lungVolume] =lungwater_ki.lungSegPythonModel(squeeze(SET(no).IM(:,:,tf,:)), [SET(no).ResolutionX SET(no).ResolutionY SET(no).SliceThickness],0);
        [mask, lungVolume] = lungwater_ki.lungSegmentation_nnUnet(SET(no).IM(:,:,tf,:), [SET(no).ResolutionX SET(no).ResolutionY SET(no).SliceThickness], tf);
        SET(no).LungWater.LungMask(:,:,1:SET(no).TSize,:)=permute(repmat(mask, [1,1, 1, SET(no).TSize]), [1 2 4 3]);
        SET(no).LungWater.LungVolume=lungVolume*ones(1,SET(no).TSize); %in ml
        lungwater_ki.lungmask2roi(no, SET(no).LungWater.LungMask,1);
        for tfs=1:SET(no).TSize
            [SET(no).LungWater.RightLungMask(:,:,tfs,:),SET(no).LungWater.LeftLungMask(:,:,tfs,:)]=lungwater_ki.rightLeftMask(no,tfs,0);
        end
    else %updates mask from ROIs
        disp('Updating manual corrections to lung mask');
        mask_new=lungwater_ki.makeLungMaskFromRoiTimeResolved(no);
        SET(no).LungWater.LungMask=mask_new;
    end

    if ~isfield(SET(no).LungWater, 'ImageTikhonov')
        disp('Normalize image');
        lambda=40.75;
        % lambda=[];
        lungwater_ki.tikhonovReg(no, lambda);
    end
    
    if ~doLiverRoi
        disp('Body mask');
        lungwater_ki.bodymask(no,bodymaskslices);

        disp('LWD calculation');
        lungwater_ki.calc_lwd(no);

    else
        disp('Liver mask');
        lungwater_ki.placeLiverROI(no);

        disp('LWD calculation');
        lungwater_ki.calc_lwd(no, doLiverRoi);
    end
    
    %regional lung water analysis
    lungwater_ki.calc_regional_lwd(no);
    lungwater_ki.postAntGradient(no, nbr_of_regions,0,filename);

    if ~isfield(SET(no).LungWater, 'TimeVector') && SET(no).TSize>1
        SET(no).LungWater.TimeVector=linspace(1/3,10, SET(no).TSize);
    elseif ~isfield(SET(no).LungWater, 'TimeVector') && SET(no).TSize==1
        SET(no).LungWater.TimeVector=0;
    end
    
     % %automatically place a background ROI in the central slice and calcualte aSNR based on a lung roi in the same slice
     % lungwater_ki.place_aSNR_ROI(no);
     % lungwater_ki.calc_aSNR(no);

     if doPlot
         tf=1;
         lungMask=squeeze(permute(SET(no).LungWater.LungMask(:,:,tf,:), [1 4 3 2]));
         lungSlice=nansum(squeeze(lungMask),3); %sum all slices to one
         lungSlice(lungSlice>0)=1; %make binary
         maskCols=find(sum(lungSlice,1)>0); %find columns with lung tissue
         sl=maskCols(round(maskCols/2)-2:round(maskCols/2)+2); %five central slices with lung tissue
         sl=1:SET(no).ZSize;
         [~,tf_peak]=max(SET(no).LungWater.DeltaLWD);

         [~, name]=fileparts(SET(no).FileName);
         str=[name ': ' SET(no).SeriesDescription];
         str_inds=strfind(str, '_');
         str(str_inds)=' ';

         h =  findobj('type','figure');
         fignbr = length(h) + 1;
         if tf_peak==1
            lungwater_ki.LWDMap_overlay(SET(no).LungWater.LWD_map(:,:,:,sl),SET(no).LungWater.LungMask(:,:,:,sl)+SET(no).LungWater.BodyMask(:,:,:,sl),SET(no).LungWater.ImageTikhonov(:,:,:,sl),str,fignbr, 0, 50, 0, tf);
         else 

            lwd_map=cat(4, SET(no).LungWater.LWD_map(:,:,tf_peak,sl), SET(no).LungWater.LWD_map(:,:,tf,sl));
            lwd_mask=cat(4, SET(no).LungWater.LungMask(:,:,tf_peak,sl), SET(no).LungWater.LungMask(:,:,tf,sl));
            lwd_norm=cat(4, SET(no).LungWater.ImageTikhonov(:,:,tf_peak,sl), SET(no).LungWater.ImageTikhonov(:,:,tf,sl));
            lungwater_ki.LWDMap_overlay(lwd_map,lwd_mask,lwd_norm,str, fignbr, 0, 50, 0, tf);

            % lungwater_ki.LWDMap_overlay(SET(no).LungWater.LWD_map(:,:,tf_peak,sl),SET(no).LungWater.LungMask(:,:,tf_peak,sl),SET(no).LungWater.ImageTikhonov(:,:,tf_peak,sl),str, fignbr, 0, 50, 0, tf);
              if SET(no).TSize>1
                 % figure(33); clf; hold on;

                 axes1 = gca; %get a handle to this axis
                 set(axes1,'OuterPosition',[0 0.5 1 0.5]);
                 axes2 = axes('Parent',figure(fignbr),'OuterPosition',[0.25 0 0.5 0.5]);

                 plot(SET(no).LungWater.TimeVector, SET(no).LungWater.LWD,'k-o'); hold on;
                 plot(SET(no).LungWater.TimeVector, SET(no).LungWater.LWD_movmean,'r-');
                 xlabel('Time (min)');
                 ylabel('Lung water density (%)');
                 legend('LWD (%)', 'Sliding mean LWD (%)')
             end
         end
     end


     disp('Post-processing finished.');
end

end


