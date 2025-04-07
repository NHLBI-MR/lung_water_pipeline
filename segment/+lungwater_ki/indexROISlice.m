function [ind] = indexROISlice(no,str,slice)
%Find index of all ROI's named str

global SET
ind=[];
for i = 1:length(SET(no).Roi)
    if isequal(SET(no).Roi(i).Name,str)
        ind =[ind i];
    end
end

if nargin >=3
    ind=ind(find([SET(no).Roi(ind).Z]==slice));
end
