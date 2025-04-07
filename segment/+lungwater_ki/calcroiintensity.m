function [m,sd,rmin,rmax,roimask] = calcroiintensity(no,roino,normalize,thisframeonly)
%Calculates intensity within a ROI (helper fcn)

global SET
if nargin < 3
    normalize = false;
end
if nargin < 4
    thisframeonly = false;
end

if thisframeonly
    tvec = SET(no).CurrentTimeFrame;
else
    tvec = SET(no).Roi(roino).T;
end

m = nan(1,SET(no).TSize);
sd = m;
rmin = m;
rmax = m;
z = SET(no).Roi(roino).Z;
for tloop=tvec
    if not(isnan(SET(no).Roi(roino).Y(1,tloop)))
        if normalize
            temp = SET(no).IM(:,:,tloop,z);
        else
            temp = calctruedata(SET(no).IM(:,:,tloop,z),no);
        end
        roimask = segment('createmask',...
            [SET(no).XSize SET(no).YSize],...
            SET(no).Roi(roino).Y(:,tloop),...
            SET(no).Roi(roino).X(:,tloop));
        ind = find(roimask);
        if ~isempty(ind)
            m(tloop) = nanmean(temp(ind));
            sd(tloop) = nanstd(temp(ind));
            rmin(tloop) = nanmin(temp(ind));
            rmax(tloop) = nanmax(temp(ind));
        end
    end
end

if thisframeonly
    m = m(tvec);
    sd = sd(tvec);
    rmin = rmin(tvec);
    rmax = rmax(tvec);
end

end

function z = calctruedata(im,no)
%-------------------------------
%Calculate true image intensities (as before Segment internal
%normalization). Uses IntensityScaling and IntensityOffset stored
%in SET structure. im is input image, and no is image stack,
%where to take the scaling from.

global SET NO

if nargin<2
    no = NO;
end

if ~isempty(SET(no).IntensityScaling)
    if not(isa(im,'int16'))
        z = im*SET(no).IntensityScaling+SET(no).IntensityOffset;
    else
        z = single(im);
    end
elseif ~isempty(SET(no).VENC) && SET(no).VENC ~= 0
    z = (im-0.5)*2*SET(no).VENC;
else
    z = im;
end
end