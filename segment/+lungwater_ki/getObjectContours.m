function [x, y]=getObjectContours(mask)

global DATA

%finds the contours around the object (only) in mask
% mask = bwareafilt(logical(mask_in), 1);
DATA.NumPoints=80;
[mi,mj] = ind2sub(size(mask),find(mask,1));

X = bwtraceboundary(mask,[mi,mj],'W');
windowWidth = 45;
polynomialOrder = 12;
polynomialOrder = 20;

if size(X,1)<windowWidth
    if mod(size(X,1),2)~=0 %odd
        windowWidth = size(X,1)-2;
    else %even
        windowWidth = size(X,1)-1;
    end
end

if windowWidth>polynomialOrder
    x = sgolayfilt([X(end-ceil(windowWidth/2):end-1,1)',X(:,1)',X(2:ceil(windowWidth/2),1)'], polynomialOrder, windowWidth);
    y = sgolayfilt([X(end-ceil(windowWidth/2):end-1,2)',X(:,2)',X(2:ceil(windowWidth/2),2)'], polynomialOrder, windowWidth);
    x = x(ceil(windowWidth/2):end-ceil(windowWidth/2)+1);
    y = y(ceil(windowWidth/2):end-ceil(windowWidth/2)+1);
    [x,y] = calcfunctions('resamplecurve',x',y',DATA.NumPoints-1);
    x=[x,x(1)]'; y=[y,y(1)]'; %first point same as last
else
    x=nan(1,DATA.NumPoints); y=nan(1,DATA.NumPoints);
end