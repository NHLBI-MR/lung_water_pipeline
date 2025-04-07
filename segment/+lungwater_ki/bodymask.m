function bodymask(no, bodymaskslices,tf)
global SET

if nargin<2
    bodymaskslices=[];
end

if nargin<3
    tf=1;
end



clear void body
mask=squeeze(SET(no).LungWater.LungMask(:,:,tf,:));

se = strel('sphere',4);
void(:,:,1,:) = imdilate(mask, se);

se = strel('sphere',8);
body(:,:,1,:) = imdilate(mask, se);
body=body-void;


%make sure there is no overlap between body and lung mask
body(body==SET(no).LungWater.LungMask(:,:,tf,:))=0;

noBodySlices=1:size(mask,3);
if isempty(bodymaskslices)
    noBodySlices(round(size(mask,3)/2)-2:round(size(mask,3)/2)+2)=[];
else
    noBodySlices(bodymaskslices)=[];
end
body(:,:,:,noBodySlices)=0*body(:,:,noBodySlices);
SET(no).LungWater.BodyMask = repmat(body,1,1,SET(no).TSize,1);
% SET(no).LungWater.BodyMask(:,:,tf,:) = body; 

bodymask2roi(no, SET(no).LungWater.BodyMask, 1);

end


function bodymask2roi(no, mask, doRegTimeRes)
%this function stores ROI's delineating the body mask
global SET DATA NO

if nargin < 3
    doRegTimeRes=0;
end

str='Body';
NO=no;

n=length(SET(no).Roi);
SET(no).RoiN=n;

DATA.Silent=true;
try
    roi('roitemplatedelete_Callback','Body');
catch
end

try
    roi('roitemplatedelete_Callback','Holes');
catch
end

try
    roi('roitemplatedelete_Callback','Void');
catch
end
DATA.Silent=false;

n=length(SET(no).Roi);

tf=1;
thistimeframeolny=1;

for slice=1:size(mask,4) %loop over slices
    cc = bwconncomp(squeeze(mask(:,:,tf,slice)),8);
    labeled = labelmatrix(cc);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [biggest,idx] =sort(numPixels,'descend');
    inds=idx(biggest>5);

    for i=1:length(inds) %loop over binary objects in mask
        m=(labeled==inds(i));
        [x, y]=lungwater_ki.getObjectContours(m);

        
        %store objects and holes as ROIs
        n = n+1;
        SET(NO).Roi(n).X = nan(length(x),SET(NO).TSize);
        SET(NO).Roi(n).Y = nan(length(y),SET(NO).TSize);
        SET(NO).Roi(n).Area = nan(1,SET(NO).TSize);
        SET(NO).Roi(n).Mean = nan(1,SET(NO).TSize);
        SET(NO).Roi(n).StD = nan(1,SET(NO).TSize);
        SET(NO).Roi(n).T = 1:SET(NO).TSize;
        if doRegTimeRes
            tfs=1:SET(NO).TSize;
        else
            tfs=1;
        end
        for tf=tfs
            SET(NO).CurrentTimeFrame=tf;

            SET(NO).Roi(n).X(:,tf)=x(:);
            SET(NO).Roi(n).Y(:,tf) = y(:);
            SET(NO).Roi(n).Z = slice;
            SET(NO).Roi(n).Sign = 1;
            % [~,area] = calcfunctions('calcroiarea',NO,n,thistimeframeolny);
            [~,area] = lungwater_ki.calcroiarea(NO,n,thistimeframeolny);
            SET(NO).Roi(n).Area(:,tf) = area;
            % [m,sd]=calcfunctions('calcroiintensity',NO,n,0,thistimeframeolny);
            [m,sd]=lungwater_ki.calcroiintensity(NO,n,0,thistimeframeolny);
            SET(NO).Roi(n).Mean(:,tf) = m;
            SET(NO).Roi(n).StD(:,tf) = sd;
            SET(NO).Roi(n).Name = str;
            SET(NO).Roi(n).LineSpec = 'r-';
            SET(NO).RoiN = n;
            SET(NO).RoiCurrent = SET(NO).RoiN;
        end
    end
end

%holes in body mask
BW5=zeros(size(mask));
for t=1:size(mask,3)
    for i=1:size(mask,4)
    BW5(:,:,t,i) = imfill(mask(:,:,t,i),'holes');
    end
end

holes=BW5-mask;

for slice=1:size(holes,4) %loop over slices
    cc = bwconncomp(squeeze(holes(:,:,tf,slice)),8);
    labeled = labelmatrix(cc);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [biggest,idx] =sort(numPixels,'descend');
    inds=idx(biggest>5);

    for i=1:length(inds) %loop over binary objects in holes
        m=(labeled==inds(i));
        [x, y]=lungwater_ki.getObjectContours(m);

        
        %store objects and holes as ROIs
        n = n+1;
        SET(NO).Roi(n).X = nan(length(x),SET(NO).TSize);
        SET(NO).Roi(n).Y = nan(length(y),SET(NO).TSize);
        SET(NO).Roi(n).Area = nan(1,SET(NO).TSize);
        SET(NO).Roi(n).Mean = nan(1,SET(NO).TSize);
        SET(NO).Roi(n).StD = nan(1,SET(NO).TSize);
        SET(NO).Roi(n).T = 1:SET(NO).TSize;
        if doRegTimeRes
            tfs=1:SET(NO).TSize;
        else
            tfs=1;
        end
        for tf=tfs
            SET(NO).CurrentTimeFrame=tf;

            SET(NO).Roi(n).X(:,tf)=x(:);
            SET(NO).Roi(n).Y(:,tf) = y(:);
            SET(NO).Roi(n).Z = slice;
            SET(NO).Roi(n).Sign = 1;
            % [~,area] = calcfunctions('calcroiarea',NO,n,thistimeframeolny);
            [~,area] = lungwater_ki.calcroiarea(NO,n,thistimeframeolny);
            SET(NO).Roi(n).Area(:,tf) = area;
            % [m,sd]=calcfunctions('calcroiintensity',NO,n,0,thistimeframeolny);
            [m,sd]=lungwater_ki.calcroiintensity(NO,n,0,thistimeframeolny);
            SET(NO).Roi(n).Mean(:,tf) = m;
            SET(NO).Roi(n).StD(:,tf) = sd;
            SET(NO).Roi(n).Name = 'Void';
            SET(NO).Roi(n).LineSpec = 'g-';
            SET(NO).RoiN = n;
            SET(NO).RoiCurrent = SET(NO).RoiN;
        end
    end
end

SET(NO).RoiCurrent = SET(NO).RoiN;

panels = find(ismember(DATA.ViewPanels,SET(no).Linked));
for p = panels
    drawfunctions('drawroi',p);
end

end