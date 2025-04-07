function calc_regional_lwd(no)
global SET


try 
    if isfield(SET(no).LungWater.Recon, 'TemporalIncr')
        dt=SET(no).LungWater.Recon.TemporalIncr;
    end
catch
    dt=1/3;
end

for tf=1:SET(no).TSize %loop over timeframes
    %pull lung and body masks
    leftLungMask=squeeze(SET(no).LungWater.LeftLungMask(:,:,tf,:));
    rightLungMask=squeeze(SET(no).LungWater.RightLungMask(:,:,tf,:));
    LWD_map=squeeze(SET(no).LungWater.LWD_map(:,:,tf,:));

    %store global LWD in struct
    SET(no).LungWater.LeftLung.LungVolume(1,tf)=sum(leftLungMask(:))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-6; %in Liter
    SET(no).LungWater.LeftLung.LWD(1,tf)=mean(LWD_map(leftLungMask==1));

    SET(no).LungWater.RightLung.LungVolume(1,tf)=sum(rightLungMask(:))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-6; %in Liter
    SET(no).LungWater.RightLung.LWD(1,tf)=mean(LWD_map(rightLungMask==1));
end

%Lung water volume
SET(no).LungWater.LeftLung.LWV=SET(no).LungWater.LeftLung.LungVolume.*SET(no).LungWater.LeftLung.LWD/100*1000; %in ml
SET(no).LungWater.RightLung.LWV=SET(no).LungWater.RightLung.LungVolume.*SET(no).LungWater.RightLung.LWD/100*1000; %in ml

%Compute and store delta LWD over time
SET(no).LungWater.LeftLung.DeltaLWD=100*(SET(no).LungWater.LeftLung.LWD-SET(no).LungWater.LeftLung.LWD(1));
SET(no).LungWater.RightLung.DeltaLWD=100*(SET(no).LungWater.RightLung.LWD-SET(no).LungWater.RightLung.LWD(1));

%Accumulation and clearance rates
%Accumulation and clearance rates
if length(SET(no).LungWater.LWD)>1
    SET(no).LungWater.LeftLung.LWD_movmean=movmean(SET(no).LungWater.LeftLung.LWD,round(1/dt));
    SET(no).LungWater.LeftLung.dLWD_dt=gradient(SET(no).LungWater.LeftLung.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.LeftLung.dLWD_dt_Max=max(SET(no).LungWater.LeftLung.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.LeftLung.dLWD_dt_Min=min(SET(no).LungWater.LeftLung.dLWD_dt); %peak clearance rate, unit %/min
    
    SET(no).LungWater.RightLung.LWD_movmean=movmean(SET(no).LungWater.RightLung.LWD,round(1/dt));
    SET(no).LungWater.RightLung.dLWD_dt=gradient(SET(no).LungWater.RightLung.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.RightLung.dLWD_dt_Max=max(SET(no).LungWater.RightLung.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.RightLung.dLWD_dt_Min=min(SET(no).LungWater.RightLung.dLWD_dt); %peak clearance rate, unit %/min
else
    SET(no).LungWater.LeftLung.LWD_movmean=[];
    SET(no).LungWater.LeftLung.dLWD_dt=[];
    SET(no).LungWater.LeftLung.dLWD_dt_Max=[];
    SET(no).LungWater.LeftLung.dLWD_dt_Min=[];
    
    SET(no).LungWater.RightLung.LWD_movmean=[];
    SET(no).LungWater.RightLung.dLWD_dt=[];
    SET(no).LungWater.RightLung.dLWD_dt_Max=[];
    SET(no).LungWater.RightLung.dLWD_dt_Min=[];
end


%Compute regional LWD (anterior, mid, and posterior)
for tf=1:SET(no).TSize %loop over timeframes

    %pull lung and body masks and permute to sagittal
    lungMask=squeeze(SET(no).LungWater.LungMask(:,:,tf,:));
    lungMask=permute(lungMask,[1,3,2]);
    lungMask(lungMask==0)=NaN;
    LWD_map=squeeze(SET(no).LungWater.LWD_map(:,:,tf,:));
    LWD_map=permute(LWD_map,[1,3,2]);


    %divide into anterior, mid , and posterior slices
    lungSlice=nansum(squeeze(lungMask),3); %sum all slices to one
    lungSlice(lungSlice>0)=1; %make binary
    maskSlices=find(sum(lungSlice,1)>0); %find slices with lung tissue
    posterior=maskSlices(1:floor(length(maskSlices)/3));
    mid=maskSlices(floor(length(maskSlices)/3)+1:2*floor(length(maskSlices)/3));
    anterior=maskSlices(2*floor(length(maskSlices)/3)+1:end);

    SET(no).LungWater.Anterior.Slices=anterior;
    SET(no).LungWater.Mid.Slices=mid;
    SET(no).LungWater.Posterior.Slices=posterior;
    SET(no).LungWater.SlicesWithLung=maskSlices;

    %regional LWD
    LWD_anterior=(LWD_map(:,anterior,:).*lungMask(:,anterior,:));
    LWD_anterior=nanmean(LWD_anterior(:));

    LWD_mid=(LWD_map(:,mid,:).*lungMask(:,mid,:));
    LWD_mid=nanmean(LWD_mid(:));

    LWD_posterior=(LWD_map(:,posterior,:).*lungMask(:,posterior,:));
    LWD_posterior=nanmean(LWD_posterior(:));

    %store regional LWD and slices in struct
    SET(no).LungWater.Anterior.LWD(1,tf)=LWD_anterior;
    SET(no).LungWater.Mid.LWD(1,tf)=LWD_mid;
    SET(no).LungWater.Posterior.LWD(1,tf)=LWD_posterior;

    SET(no).LungWater.Anterior.LWV(1,tf)=LWD_anterior/100*sum(sum(sum(lungMask(:,anterior,:)==1)))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-3; %in ml
    SET(no).LungWater.Mid.LWV(1,tf)=(LWD_mid/100)*sum(sum(sum(lungMask(:,mid,:)==1)))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-3; %in ml
    SET(no).LungWater.Posterior.LWV(1,tf)=LWD_posterior/100*sum(sum(sum(lungMask(:,posterior,:)==1)))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-3; %in ml

end

SET(no).LungWater.Anterior.DeltaLWD=100*(SET(no).LungWater.Anterior.LWD-SET(no).LungWater.Anterior.LWD(1));
SET(no).LungWater.Mid.DeltaLWD=100*(SET(no).LungWater.Mid.LWD-SET(no).LungWater.Mid.LWD(1));
SET(no).LungWater.Posterior.DeltaLWD=100*(SET(no).LungWater.Posterior.LWD-SET(no).LungWater.Posterior.LWD(1));

SET(no).LungWater.AnteriorPosteriorLength=(length(SET(no).LungWater.SlicesWithLung).*SET(no).SliceThickness./10); %cm
SET(no).LungWater.AnteriorPosteriorGradient=(SET(no).LungWater.Posterior.LWD(1)-SET(no).LungWater.Anterior.LWD(1))/SET(no).LungWater.AnteriorPosteriorLength; %unit %/cm

%regional accumulation and clearance rates
if length(SET(no).LungWater.LWD)>1
    SET(no).LungWater.Anterior.LWD_movmean=movmean(SET(no).LungWater.Anterior.LWD,round(1/dt));
    SET(no).LungWater.Anterior.dLWD_dt=gradient(SET(no).LungWater.Anterior.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.Anterior.dLWD_dt_Max=max(SET(no).LungWater.Anterior.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.Anterior.dLWD_dt_Min=min(SET(no).LungWater.Anterior.dLWD_dt); %peak clearance rate, unit %/min

    SET(no).LungWater.Mid.LWD_movmean=movmean(SET(no).LungWater.Mid.LWD,round(1/dt));
    SET(no).LungWater.Mid.dLWD_dt=gradient(SET(no).LungWater.Mid.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.Mid.dLWD_dt_Max=max(SET(no).LungWater.Mid.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.Mid.dLWD_dt_Min=min(SET(no).LungWater.Mid.dLWD_dt); %peak clearance rate, unit %/min

    SET(no).LungWater.Posterior.LWD_movmean=movmean(SET(no).LungWater.Posterior.LWD,round(1/dt));
    SET(no).LungWater.Posterior.dLWD_dt=gradient(SET(no).LungWater.Posterior.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.Posterior.dLWD_dt_Max=max(SET(no).LungWater.Posterior.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.Posterior.dLWD_dt_Min=min(SET(no).LungWater.Posterior.dLWD_dt); %peak clearance rate, unit %/min
else
    SET(no).LungWater.Anterior.LWD_movmean=[];
    SET(no).LungWater.Anterior.dLWD_dt=[];
    SET(no).LungWater.Anterior.dLWD_dt_Max=[];
    SET(no).LungWater.Anterior.dLWD_dt_Min=[];

    SET(no).LungWater.Mid.LWD_movmean=[];
    SET(no).LungWater.Mid.dLWD_dt=[];
    SET(no).LungWater.Mid.dLWD_dt_Max=[];
    SET(no).LungWater.Mid.dLWD_dt_Min=[];

    SET(no).LungWater.Posterior.LWD_movmean=[];
    SET(no).LungWater.Posterior.dLWD_dt=[];
    SET(no).LungWater.Posterior.dLWD_dt_Max=[];
    SET(no).LungWater.Posterior.dLWD_dt_Min=[];
end

%Compute regional LWD (apical, central, basal)
for tf=1:SET(no).TSize %loop over timeframes

    %pull lung and body masks and permute to sagittal
    lungMask=squeeze(SET(no).LungWater.LungMask(:,:,tf,:));
    lungMask(lungMask==0)=NaN;
    LWD_map=squeeze(SET(no).LungWater.LWD_map(:,:,tf,:));

    %divide into apical, central , and basal slices
    lungSlice=nansum(squeeze(lungMask),3); %sum all slices to one
    lungSlice(lungSlice>0)=1; %make binary
    maskRows=find(sum(lungSlice,2)>0); %find rows with lung tissue
    apical=maskRows(1:floor(length(maskRows)/3));
    central=maskRows(floor(length(maskRows)/3)+1:2*floor(length(maskRows)/3));
    basal=maskRows(2*floor(length(maskRows)/3)+1:end);

    SET(no).LungWater.Apical.Rows=apical;
    SET(no).LungWater.Central.Rows=central;
    SET(no).LungWater.Basal.Rows=basal;

    %regional LWD
    LWD_apical=(LWD_map(:,apical,:).*lungMask(:,apical,:));
    LWD_apical=nanmean(LWD_apical(:));

    LWD_central=(LWD_map(:,central,:).*lungMask(:,central,:));
    LWD_central=nanmean(LWD_central(:));

    LWD_basal=(LWD_map(:,basal,:).*lungMask(:,basal,:));
    LWD_basal=nanmean(LWD_basal(:));

    %store regional LWD and slices in struct
    SET(no).LungWater.Apical.LWD(1,tf)=LWD_apical;
    SET(no).LungWater.Central.LWD(1,tf)=LWD_central;
    SET(no).LungWater.Basal.LWD(1,tf)=LWD_basal;

    SET(no).LungWater.Apical.LWV(1,tf)=LWD_apical/100*sum(sum(sum(lungMask(:,apical,:)==1)))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-3; %in ml
    SET(no).LungWater.Central.LWV(1,tf)=(LWD_central/100)*sum(sum(sum(lungMask(:,central,:)==1)))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-3; %in ml
    SET(no).LungWater.Basal.LWV(1,tf)=LWD_basal/100*sum(sum(sum(lungMask(:,basal,:)==1)))*(SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap))*10^-3; %in ml

end

SET(no).LungWater.Apical.DeltaLWD=100*(SET(no).LungWater.Apical.LWD-SET(no).LungWater.Apical.LWD(1));
SET(no).LungWater.Central.DeltaLWD=100*(SET(no).LungWater.Central.LWD-SET(no).LungWater.Central.LWD(1));
SET(no).LungWater.Basal.DeltaLWD=100*(SET(no).LungWater.Basal.LWD-SET(no).LungWater.Basal.LWD(1));

%regional accumulation and clearance rates
if length(SET(no).LungWater.LWD)>1
    SET(no).LungWater.Apical.LWD_movmean=movmean(SET(no).LungWater.Apical.LWD,round(1/dt));
    SET(no).LungWater.Apical.dLWD_dt=gradient(SET(no).LungWater.Apical.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.Apical.dLWD_dt_Max=max(SET(no).LungWater.Apical.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.Apical.dLWD_dt_Min=min(SET(no).LungWater.Apical.dLWD_dt); %peak clearance rate, unit %/min

    SET(no).LungWater.Central.LWD_movmean=movmean(SET(no).LungWater.Central.LWD,round(1/dt));
    SET(no).LungWater.Central.dLWD_dt=gradient(SET(no).LungWater.Central.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.Central.dLWD_dt_Max=max(SET(no).LungWater.Central.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.Central.dLWD_dt_Min=min(SET(no).LungWater.Central.dLWD_dt); %peak clearance rate, unit %/min

    SET(no).LungWater.Basal.LWD_movmean=movmean(SET(no).LungWater.Basal.LWD,round(1/dt));
    SET(no).LungWater.Basal.dLWD_dt=gradient(SET(no).LungWater.Basal.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.Basal.dLWD_dt_Max=max(SET(no).LungWater.Basal.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.Basal.dLWD_dt_Min=min(SET(no).LungWater.Basal.dLWD_dt); %peak clearance rate, unit %/min
else
    SET(no).LungWater.Apical.LWD_movmean=[];
    SET(no).LungWater.Apical.dLWD_dt=[];
    SET(no).LungWater.Apical.dLWD_dt_Max=[];
    SET(no).LungWater.Apical.dLWD_dt_Min=[];

    SET(no).LungWater.Central.LWD_movmean=[];
    SET(no).LungWater.Central.dLWD_dt=[];
    SET(no).LungWater.Central.dLWD_dt_Max=[];
    SET(no).LungWater.Central.dLWD_dt_Min=[];

    SET(no).LungWater.Basal.LWD_movmean=[];
    SET(no).LungWater.Basal.dLWD_dt=[];
    SET(no).LungWater.Basal.dLWD_dt_Max=[];
    SET(no).LungWater.Basal.dLWD_dt_Min=[];
end
