function lungmask2roi(no, mask, doRegTimeRes)
global SET DATA NO

if nargin < 3
    doRegTimeRes=0;
end

str='Lung';
NO=no;

DATA.Silent=true;
roi('roitemplatedelete_Callback',str);
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
%         if doRegTimeRes
            tfs=1:SET(NO).TSize;
%         else
%             tfs=1;
%         end
        for tf=tfs
            SET(NO).CurrentTimeFrame=tf;

            SET(NO).Roi(n).X(:,tf)=x(:);
            SET(NO).Roi(n).Y(:,tf) = y(:);
            %             SET(NO).Roi(n).T(:,tf) = tf;
            SET(NO).Roi(n).Z = slice;
            SET(NO).Roi(n).Sign = 1;
            [~,area] = calcfunctions('calcroiarea',NO,n,thistimeframeolny);
            SET(NO).Roi(n).Area(:,tf) = area;
            [m,sd]=calcfunctions('calcroiintensity',NO,n,0,thistimeframeolny);
            SET(NO).Roi(n).Mean(:,tf) = m;
            SET(NO).Roi(n).StD(:,tf) = sd;
            SET(NO).Roi(n).Name = str;
            SET(NO).Roi(n).LineSpec = 'b-';
            SET(NO).RoiN = n;
            SET(NO).RoiCurrent = SET(NO).RoiN;
        end
    end
     tf=1;
end



SET(NO).RoiCurrent = SET(NO).RoiN;

for rois=lungwater_ki.indexROISlice(no,'Lung')
    SET(NO).RoiCurrent = rois;
    roi('expandcontract_Callback',2);
end

panels = find(ismember(DATA.ViewPanels,SET(no).Linked));
for p = panels
    drawfunctions('drawroi',p);
end

