function calc_lwd(no, doLiverRoi)

if nargin <2
  doLiverRoi=0;  
end
%ADD regional (right and left lung, anterior, mid, posterior)

global SET

%pull normalized image
image=SET(no).LungWater.ImageTikhonov;
PdScaling=70; %assume 70% musculosceletal water density

try 
    if isfield(SET(no).LungWater.Recon, 'TemporalIncr')
        dt=SET(no).LungWater.Recon.TemporalIncr;
    end
catch
    dt=1/3;
end


for tf=1:SET(no).TSize %loop over timeframes
    %make image 3D
    im=squeeze(double(image(:,:,tf,:)));

    %pull lung and body masks
    lungMask=squeeze(SET(no).LungWater.LungMask(:,:,tf,:));
    if ~doLiverRoi
        bodyMask=squeeze(SET(no).LungWater.BodyMask(:,:,tf,:));
    else
        bodyMask=squeeze(lungwater_ki.makeLungMaskFromRoiTimeResolved(no,'Liver'));
    end

    %lung mask
    lungMask(lungMask==0)=NaN; %make zero elements in lung mask to NaN for easier computation
    lungSI=nanmean(im(lungMask==1)); %mean lung signal intensity
    bodySI=median(im(bodyMask==1)); %median body signal (median, not mean, to avoid trachea and hyper ehnaced SIs)

    %global LWD
    LWD_map=PdScaling*(im)/bodySI; %pixel-wise LWD calculation
    LWD=nanmean(LWD_map(:).*lungMask(:)); %global LWD

    %store global LWD in struct
%     SET(no).LungWater.LungSI(1,tf)=lungSI;
%     SET(no).LungWater.BodySI(1,tf)=bodySI;
    SET(no).LungWater.LWD_map(:,:,tf,:)=LWD_map;
    SET(no).LungWater.LWD(1,tf)=LWD;
end

%Lung water volume
SET(no).LungWater.LWV=SET(no).LungWater.LungVolume.*SET(no).LungWater.LWD/100*1000; %in ml

%Compute and store delta LWD over time
SET(no).LungWater.DeltaLWD=100*(SET(no).LungWater.LWD-SET(no).LungWater.LWD(1)); % absolute difference, not in percent of baseline value

%Accumulation and clearance rates
if length(SET(no).LungWater.LWD)>1
    SET(no).LungWater.LWD_movmean=movmean(SET(no).LungWater.LWD,round(1/dt));
    SET(no).LungWater.dLWD_dt=gradient(SET(no).LungWater.LWD_movmean,dt); %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.dLWD_dt_Max=max(SET(no).LungWater.dLWD_dt); %peak accumulation rate, unit %/min
    SET(no).LungWater.dLWD_dt_Min=min(SET(no).LungWater.dLWD_dt); %peak clearance rate, unit %/min
else
    SET(no).LungWater.LWD_movmean=[];
    SET(no).LungWater.dLWD_dt=[]; %derivative of 1-minute sliding mean. dt=1/3 min = 20 s so 1 min = 3 timeframs
    SET(no).LungWater.dLWD_dt_Max=[]; %peak accumulation rate, unit %/min
    SET(no).LungWater.dLWD_dt_Min=[]; %peak clearance rate, unit %/min
end
