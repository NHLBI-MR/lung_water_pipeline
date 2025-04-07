function place_aSNR_ROI(nos)

%add check if ROI overlaps with lung segmentation, in that case move ROI
%down. I antecipate this function will be phased out.


global SET NO DATA

for no=nos
NO=no;
roi('roitemplatedelete_Callback','aSNR'); %delete any Roi named aSNR

tSize=SET(no).TSize;

midSlice=round(SET(no).ZSize/2);    
col=SET(no).YSize-10; %column for center of background ROI
row=7; %find row for center of background ROI
    
n = size(SET(NO).Roi,2)+1;
SET(NO).Roi(n).X = nan(80,tSize);
SET(NO).Roi(n).Y = nan(80,tSize);
SET(NO).Roi(n).Z = midSlice;
SET(NO).Roi(n).Area = nan(1,tSize);
SET(NO).Roi(n).Mean = nan(1,tSize);
SET(NO).Roi(n).StD = nan(1,tSize);
SET(NO).Roi(n).T = 1:tSize;

for tf=1:tSize %loop over timeframes

    %coordinates for center of ROI
    r0=10; %radius for liver ROI
    np = DATA.NumPoints-1;
    omega = ((2*pi/np)*(1:np))';
    rx = r0/SET(no).ResolutionX;
    ry = r0/SET(no).ResolutionY;
    X=row;%10;
    Y=col;
    
    %create circular mask
    x = repmat(rx*sin(omega)+X,[1 1 1]);
    y = repmat(ry*cos(omega)+Y,[1 1 1]);
    x = [x ; x(1,:)];
    y = [y ; y(1,:)];
    
    
    %store mask as ROI
    NO=no;
    SET(NO).CurrentTimeFrame=tf;
    
%     SET(NO).Roi(n).X = nan(length(x),tSize);
%     SET(NO).Roi(n).Y = nan(length(y),tSize);
%     SET(NO).Roi(n).Area = nan(1,tSize);
%     SET(NO).Roi(n).Mean = nan(1,tSize);
%     SET(NO).Roi(n).StD = nan(1,tSize);
    
    SET(NO).Roi(n).X(:,tf) = x(:);
    SET(NO).Roi(n).Y(:,tf) = y(:);
%     SET(NO).Roi(n).T = 1:tf;
    SET(NO).Roi(n).Z = midSlice;
    SET(NO).Roi(n).Sign = 1;
    % [~,area] = calcfunctions('calcroiarea',NO,n);
    [~,area] = lungwater_ki.calcroiarea(NO,n);
    SET(NO).Roi(n).Area = area;
    % [m,sd]=calcfunctions('calcroiintensity',NO,n);
    [m,sd]=lungwater_ki.calcroiintensity(NO,n);
    SET(NO).Roi(n).Mean = m;
    SET(NO).Roi(n).StD = sd;
    SET(NO).Roi(n).Name = 'aSNR';
    SET(NO).Roi(n).LineSpec = 'r-';
    SET(NO).RoiN = n;
    SET(NO).RoiCurrent = SET(NO).RoiN;
end

%store in SET struct
SET(no).LungWater.aSNRMask=[]; 
SET(no).LungWater.aSNRMask=lungwater_ki.makeLungMaskFromRoiTimeResolved(no,'aSNR');

SET(NO).RoiCurrent = SET(NO).RoiN;

panels = find(ismember(DATA.ViewPanels,SET(no).Linked));
for p = panels
    drawfunctions('drawroi',p);
end

end