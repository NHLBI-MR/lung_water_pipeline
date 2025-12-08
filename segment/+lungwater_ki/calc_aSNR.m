function calc_aSNR(nos)

global SET NO
%find aSNR background ROI

for no=nos
    % NO=no;
    indBackground = lungwater_ki.indexROISlice(no,'aSNR');

    slice=SET(no).Roi(indBackground).Z;

    indLung = lungwater_ki.indexROISlice(no,'Lung', slice);
    aSNR=(SET(no).Roi(indLung(end)).Mean)./(SET(no).Roi(indBackground).StD);

    SET(no).LungWater.aSNR.aSNR=aSNR;
    SET(no).LungWater.aSNR.LungSIMean=SET(no).Roi(indLung(end)).Mean;
    SET(no).LungWater.aSNR.BackgroundSIStD=SET(no).Roi(indBackground).StD;

end

end