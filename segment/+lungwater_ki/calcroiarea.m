function [meanarea,area]=calcroiarea(no,roino,thisframeonly) %
%--------------------------
%Calculates roi area (helper fcn)

global SET

if nargin < 3
  thisframeonly = false;
end
if thisframeonly
  tvec = SET(no).CurrentTimeFrame;
else
  tvec = SET(no).Roi(roino).T;
end

area=nan(1,SET(no).TSize);
for tloop=tvec
  if not(isnan(SET(no).Roi(roino).Y(1,tloop)))
    area(tloop)= (1/100)*stablepolyarea(...
      SET(no).ResolutionY*SET(no).Roi(roino).Y(:,tloop),...
      SET(no).ResolutionX*SET(no).Roi(roino).X(:,tloop));
  end
end
if thisframeonly
  area = area(tvec);
end
meanarea=mynanmean(area);