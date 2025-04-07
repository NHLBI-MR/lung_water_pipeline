function out=exportDynamicLWDSegment(nos, doHeader)

global SET

output=cell(size(nos,1)+1,80);

tfs=([SET(nos).TSize]);
row=1;

for loop=nos
    no=nos(loop);

    output{row,1}=SET(1).PatientInfo.AcquisitionDate; %date
    [path,name] = fileparts(SET(1).FileName);
    output{row,2}=path;
    output{row,3}=name;
   
    %output{row,2}=SET(1).FileName; %subject ID
    output{row,8}=SET(1).PatientInfo.Sex; 
    output{row,9}=SET(1).PatientInfo.Age; 
    output{row,10}=SET(1).PatientInfo.Length; 
    output{row,11}=SET(1).PatientInfo.Weight; 
    output{row,12}=SET(1).PatientInfo.BSA; 
    output{row,13}=SET(no).SeriesDescription; 
    output{row,14}=no; %repetition nbr
    output{row,15}=lungwater_ki.segmenttime2clock(no); %acquisiton time

    n = 15;
        
    output{row,n+2}=SET(no).LungWater.LWD(1);
    output{row,n+3}=max(SET(no).LungWater.LWD);
    output{row,n+4}=SET(no).LungWater.LWD(end);
    output{row,n+5}=max(SET(no).LungWater.DeltaLWD);
    output{row,n+6}=SET(no).LungWater.dLWD_dt_Max; 
    output{row,n+7}=SET(no).LungWater.dLWD_dt_Min; 

    output{row,n+16}=SET(no).LungWater.RightLung.LWD(1);
    output{row,n+17}=max(SET(no).LungWater.RightLung.LWD);
    output{row,n+18}=SET(no).LungWater.RightLung.LWD(end);
    output{row,n+19}=max(SET(no).LungWater.RightLung.DeltaLWD);
    output{row,n+20}=SET(no).LungWater.RightLung.dLWD_dt_Max; 
    output{row,n+21}=SET(no).LungWater.RightLung.dLWD_dt_Min; 

    output{row,n+23}=SET(no).LungWater.LeftLung.LWD(1);
    output{row,n+24}=max(SET(no).LungWater.LeftLung.LWD);
    output{row,n+25}=SET(no).LungWater.LeftLung.LWD(end);
    output{row,n+26}=max(SET(no).LungWater.LeftLung.DeltaLWD);
    output{row,n+27}=SET(no).LungWater.LeftLung.dLWD_dt_Max; 
    output{row,n+28}=SET(no).LungWater.LeftLung.dLWD_dt_Min; 

    output{row,n+30}=SET(no).LungWater.Anterior.LWD(1);
    output{row,n+31}=max(SET(no).LungWater.Anterior.LWD);
    output{row,n+32}=SET(no).LungWater.Anterior.LWD(end);
    output{row,n+33}=max(SET(no).LungWater.Anterior.DeltaLWD);
    output{row,n+34}=SET(no).LungWater.Anterior.dLWD_dt_Max; 
    output{row,n+35}=SET(no).LungWater.Anterior.dLWD_dt_Min; 

    output{row,n+37}=SET(no).LungWater.Mid.LWD(1);
    output{row,n+38}=max(SET(no).LungWater.Mid.LWD);
    output{row,n+39}=SET(no).LungWater.Mid.LWD(end);
    output{row,n+40}=max(SET(no).LungWater.Mid.DeltaLWD);
    output{row,n+41}=SET(no).LungWater.Mid.dLWD_dt_Max; 
    output{row,n+42}=SET(no).LungWater.Mid.dLWD_dt_Min; 
    
    output{row,n+44}=SET(no).LungWater.Posterior.LWD(1);
    output{row,n+45}=max(SET(no).LungWater.Posterior.LWD);
    output{row,n+46}=SET(no).LungWater.Posterior.LWD(end);
    output{row,n+47}=max(SET(no).LungWater.Posterior.DeltaLWD);
    output{row,n+48}=SET(no).LungWater.Posterior.dLWD_dt_Max; 
    output{row,n+49}=SET(no).LungWater.Posterior.dLWD_dt_Min; 

    output{row,n+51}=SET(no).LungWater.Apical.LWD(1);
    output{row,n+52}=max(SET(no).LungWater.Apical.LWD);
    output{row,n+53}=SET(no).LungWater.Apical.LWD(end);
    output{row,n+54}=max(SET(no).LungWater.Apical.DeltaLWD);
    output{row,n+55}=SET(no).LungWater.Apical.dLWD_dt_Max; 
    output{row,n+56}=SET(no).LungWater.Apical.dLWD_dt_Min; 

    output{row,n+58}=SET(no).LungWater.Central.LWD(1);
    output{row,n+59}=max(SET(no).LungWater.Central.LWD);
    output{row,n+60}=SET(no).LungWater.Central.LWD(end);
    output{row,n+61}=max(SET(no).LungWater.Central.DeltaLWD);
    output{row,n+62}=SET(no).LungWater.Central.dLWD_dt_Max; 
    output{row,n+63}=SET(no).LungWater.Central.dLWD_dt_Min; 
    
    output{row,n+65}=SET(no).LungWater.Basal.LWD(1);
    output{row,n+66}=max(SET(no).LungWater.Basal.LWD);
    output{row,n+67}=SET(no).LungWater.Basal.LWD(end);
    output{row,n+68}=max(SET(no).LungWater.Basal.DeltaLWD);
    output{row,n+69}=SET(no).LungWater.Basal.dLWD_dt_Max; 
    output{row,n+70}=SET(no).LungWater.Basal.dLWD_dt_Min; 

    [~, maxInd]=max(SET(no).LungWater.LWD); %time point for peak global LWD
    output{row,n+72}=maxInd;
    output{row,n+73}=SET(no).LungWater.LWV(maxInd);
    output{row,n+74}=SET(no).LungWater.Anterior.LWD(maxInd);
    output{row,n+75}=SET(no).LungWater.Mid.LWD(maxInd);
    output{row,n+76}=SET(no).LungWater.Posterior.LWD(maxInd); 
    output{row,n+77}=SET(no).LungWater.Apical.LWD(maxInd);
    output{row,n+78}=SET(no).LungWater.Central.LWD(maxInd);
    output{row,n+79}=SET(no).LungWater.Basal.LWD(maxInd); 

    output{row,n+81}=SET(no).LungWater.Posterior.LWD(1)-SET(no).LungWater.Anterior.LWD(1);
    output{row,n+82}=SET(no).LungWater.Posterior.LWD(maxInd)-SET(no).LungWater.Anterior.LWD(maxInd);
    output{row,n+83}=SET(no).LungWater.AnteriorPosteriorLength;
    output{row,n+84}=SET(no).LungWater.AnteriorPosteriorGradient;
    output{row,n+85}=(SET(no).LungWater.Posterior.LWD(maxInd)-SET(no).LungWater.Anterior.LWD(maxInd))/SET(no).LungWater.AnteriorPosteriorLength; %unit %/cm

    for tf=1:tfs(loop)
        output{row+tf-1,n+9}=SET(no).LungWater.TimeVector(tf)+SET(no).LungWater.Recon.TemporalFootprint; 
        output{row+tf-1,n+10}=SET(no).LungWater.LungVolume(tf);
        output{row+tf-1,n+11}=SET(no).LungWater.LWD(tf);
        output{row+tf-1,n+12}=SET(no).LungWater.LWV(tf); 
        if length(SET(no).LungWater.LWD)>1
            output{row+tf-1,n+13}=SET(no).LungWater.dLWD_dt(tf); 
            output{row+tf-1,n+14}=SET(no).LungWater.DeltaLWD(tf);
            output{row+tf-1,n+15}=SET(no).LungWater.LWD_movmean(tf);
        end
        % output{row+tf-1,n+88}=SET(no).LungWater.aSNR.aSNR(tf);
    end

    output{row,n+90}=SET(no).LungWater.sliceProfilePosteriorAnterior(1);
    output{row,n+91}=SET(no).LungWater.sliceProfilePosteriorAnterior(2);
    output{row,n+92}=SET(no).LungWater.sliceProfilePosteriorAnterior(3);
    output{row,n+93}=SET(no).LungWater.sliceProfilePosteriorAnterior(4);
    output{row,n+94}=SET(no).LungWater.sliceProfilePosteriorAnterior(5);
    output{row,n+95}=SET(no).LungWater.sliceProfilePosteriorAnterior(6);
    output{row,n+96}=SET(no).LungWater.sliceProfilePosteriorAnterior(7);
    output{row,n+97}=SET(no).LungWater.sliceProfilePosteriorAnterior(8);
    output{row,n+98}=SET(no).LungWater.sliceProfilePosteriorAnterior(9);
    output{row,n+99}=SET(no).LungWater.sliceProfilePosteriorAnterior(10);

    row=row+tfs(loop)+1;
end
output{row,1}='';

if doHeader

    header=cell(1,size(output,2));
    
    header{1,1}='Acquisition Date';
    header{1,2}='Path';
    header{1,3}='Subject ID';
    header{1,4}='Facility ID';
    header{1,5}='Patient/Healthy volunteer';
    header{1,6}='Comment';

    header{1,8}='Sex';
    header{1,9}='Age (years)';
    header{1,10}='Length (cm)';
    header{1,11}='Weight (kg)';
    header{1,12}='BSA (m^2)';
    header{1,13}='Protocol';
    header{1,14}='Image stack no';
    header{1,15}='Acquistion time';

    header{1,n+2}='Baseline LWD (%)';
    header{1,n+3}='Peak LWD (%)';
    header{1,n+4}='End LWD (%)';
    header{1,n+5}='Peak delta LWD (%)';
    header{1,n+6}='dLWD/dt max (%/min)';
    header{1,n+7}='dLWD/dt min (%/min)';

    header{1,n+9}='Time (min)';
    header{1,n+10}='Lung Volume (L)';
    header{1,n+11}='LWD (%)';
    header{1,n+12}='LWV (ml)';
    header{1,n+13}='dLWD/dt (%/min)';
    header{1,n+14}='Delta LWD (%)';
    header{1,n+15}='Sliding mean LWD (%)';

    header{1,n+16}='Right Lung Baseline LWD (%)';
    header{1,n+17}='Right Lung Peak LWD (%)';
    header{1,n+18}='Right Lung End LWD (%)';
    header{1,n+19}='Right Lung Peak delta LWD (%)';
    header{1,n+20}='Right Lung dLWD/dt max (%/min)';
    header{1,n+21}='Right Lung dLWD/dt min (%/min)';

    header{1,n+23}='Left Lung Baseline LWD (%)';
    header{1,n+24}='Left Lung Peak LWD (%)';
    header{1,n+25}='Left Lung End LWD (%)';
    header{1,n+26}='Left Lung Peak delta LWD (%)';
    header{1,n+27}='Left Lung dLWD/dt max (%/min)';
    header{1,n+28}='Left Lung dLWD/dt min (%/min)';

    header{1,n+30}='Anterior Baseline LWD (%)';
    header{1,n+31}='Anterior Peak LWD (%)';
    header{1,n+32}='Anterior End LWD (%)';
    header{1,n+33}='Anterior Peak delta LWD (%)';
    header{1,n+34}='Anterior dLWD/dt max (%/min)';
    header{1,n+35}='Anterior dLWD/dt min (%/min)';

    header{1,n+37}='Mid Baseline LWD (%)';
    header{1,n+38}='Mid Peak LWD (%)';
    header{1,n+39}='Mid End LWD (%)';
    header{1,n+40}='Mid Peak delta LWD (%)';
    header{1,n+41}='Mid dLWD/dt max (%/min)';
    header{1,n+42}='Mid dLWD/dt min (%/min)';

    header{1,n+44}='Posterior Baseline LWD (%)';
    header{1,n+45}='Posterior Peak LWD (%)';
    header{1,n+46}='Posterior End LWD (%)';
    header{1,n+47}='Posterior Peak delta LWD (%)';
    header{1,n+48}='Posterior dLWD/dt max (%/min)';
    header{1,n+49}='Posterior dLWD/dt min (%/min)';

    header{1,n+51}='Apical Baseline LWD (%)';
    header{1,n+52}='Apical Peak LWD (%)';
    header{1,n+53}='Apical End LWD (%)';
    header{1,n+54}='Apical Peak delta LWD (%)';
    header{1,n+55}='Apical dLWD/dt max (%/min)';
    header{1,n+56}='Apical dLWD/dt min (%/min)';

    header{1,n+58}='Center Baseline LWD (%)';
    header{1,n+59}='Center Peak LWD (%)';
    header{1,n+60}='Center End LWD (%)';
    header{1,n+61}='Center Peak delta LWD (%)';
    header{1,n+62}='Center dLWD/dt max (%/min)';
    header{1,n+63}='Center dLWD/dt min (%/min)';

    header{1,n+65}='Basal Baseline LWD (%)';
    header{1,n+66}='Basal Peak LWD (%)';
    header{1,n+67}='Basal End LWD (%)';
    header{1,n+68}='Basal Peak delta LWD (%)';
    header{1,n+69}='Basal dLWD/dt max (%/min)';
    header{1,n+70}='Basal dLWD/dt min (%/min)';


    header{1,n+72}='Timeframe of peak delta LWD';
    header{1,n+73}='LWV at peak delta LWD (ml)';
    header{1,n+74}='Anterior LWD at peak delta LWD (%)';
    header{1,n+75}='Mid LWD at peak delta LWD (%)';
    header{1,n+76}='Posterior LWD at peak delta LWD (%)';
    header{1,n+77}='Apical LWD at peak delta LWD (%)';
    header{1,n+78}='Central LWD at peak delta LWD (%)';
    header{1,n+79}='Basal LWD at peak delta LWD (%)';


    header{1,n+81}='Posterior-Anterior Baseline LWD (%)';
    header{1,n+82}='Posterior-Anterior at peak delta LWD (%)';
    header{1,n+83}='Anteroposterior lung width (cm)';
    header{1,n+84}='Baseline Anterior-Posterior Distribution Gradient (%/cm)';
    header{1,n+85}='Peak LWD Anterior-Posterior Distribution Gradient (%/cm)';

    % header{1,n+88}='aSNR';

    header{1,n+90}='LWD slice profile bin 1 (%)';
    header{1,n+91}='LWD slice profile bin 2 (%)';
    header{1,n+92}='LWD slice profile bin 3 (%)';
    header{1,n+93}='LWD slice profile bin 4 (%)';
    header{1,n+94}='LWD slice profile bin 5 (%)';
    header{1,n+95}='LWD slice profile bin 6 (%)';
    header{1,n+96}='LWD slice profile bin 7 (%)';
    header{1,n+97}='LWD slice profile bin 8 (%)';
    header{1,n+98}='LWD slice profile bin 9 (%)';
    header{1,n+99}='LWD slice profile bin 10 (%)';

    out=[header;output];

else
    out=output;
end