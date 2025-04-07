function placeLiverROI(no, prone)

%add check if ROI overlaps with lung segmentation, in that case move ROI
%down. I antecipate this function will be phased out.

if nargin<2
   prone=0;
end

global SET NO DATA
NO=no;

roi('roitemplatedelete_Callback','Liver'); %delete any Roi named Liver

tSize=SET(no).TSize;

for tf=1:tSize %loop over timeframes
    if tf==1
        lungMask=squeeze(SET(no).LungWater.LungMask(:,:,tf,:));
        
        dilatelevel=2;
        for sl=1:size(lungMask,3)
            lungMask(:,:,sl) = imerode(lungMask(:,:,sl),strel('disk',dilatelevel));
        end

        rightLungMask=squeeze(SET(no).LungWater.RightLungMask);
        
        if prone %swap things arounf if image is taken in prone position
            rightLungMask=squeeze(SET(no).LungWater.LeftLungMask);
            shiftConstant=2/3;
        else
            shiftConstant=1/3;
        end
        
%         if sum(rightLungMask(:))<=100
%             rightLungMask=lungMask;
%         end
%         
        %find mid lung slice using centroid
        cc = bwconncomp(rightLungMask,26);
        S = regionprops(cc,'Centroid');
        mid=round(S(1).Centroid);
        midSlice=mid(3);
        
        %find middle column of 3D right lung mask
        cols=find(sum(sum(rightLungMask(:,:,midSlice),3))>0);
        midCol=cols(round(shiftConstant*length(cols)));

        %find bottom middle lung
        bottomMidRow=find(rightLungMask(:,midCol,midSlice)>0, 1, 'last'); %max(find(rightLungMask(:,midCol,midSlice)>0));
        
    end
    
    %coordinates for center of ROI
    r0=20; %radius for liver ROI
    n = DATA.NumPoints-1;
    omega = ((2*pi/n)*(1:n))';
    rx = r0/SET(no).ResolutionX;
    ry = r0/SET(no).ResolutionY;
    X=bottomMidRow+ry+10;%10;
    Y=midCol;
    
    %create circular mask
    x = repmat(rx*sin(omega)+X,[1 1 1]);
    y = repmat(ry*cos(omega)+Y,[1 1 1]);
    x = [x ; x(1,:)];
    y = [y ; y(1,:)];
    
    
    %store mask as ROI
    NO=no;
    SET(NO).CurrentTimeFrame=tf;
    n = size(SET(NO).Roi,2)+1;
    SET(NO).Roi(n).X = nan(length(x),tSize);
    SET(NO).Roi(n).Y = nan(length(y),tSize);
    SET(NO).Roi(n).Area = nan(1,tSize);
    SET(NO).Roi(n).Mean = nan(1,tSize);
    SET(NO).Roi(n).StD = nan(1,tSize);
    
    SET(NO).Roi(n).X(:,tf) = x(:);
    SET(NO).Roi(n).Y(:,tf) = y(:);
    SET(NO).Roi(n).T = tf;
    SET(NO).Roi(n).Z = midSlice;
    SET(NO).Roi(n).Sign = 1;
    % [~,area] = calcfunctions('calcroiarea',NO,n);
    [~,area] = lungwater_ki.calcroiarea(NO,n);
    SET(NO).Roi(n).Area = area;
    % [m,sd]=calcfunctions('calcroiintensity',NO,n);
    [m,sd]=lungwater_ki.calcroiintensity(NO,n);
    SET(NO).Roi(n).Mean = m;
    SET(NO).Roi(n).StD = sd;
    SET(NO).Roi(n).Name = 'Liver';
    SET(NO).Roi(n).LineSpec = 'r-';
    SET(NO).RoiN = n;
    SET(NO).RoiCurrent = SET(NO).RoiN;
end

%store in SET struct
SET(no).LungWater.LiverMask=[]; 
SET(no).LungWater.LiverMask=lungs.makeLungMaskFromRoiTimeResolved(no,'Liver');
