function clockstr = segmenttime2clock(no)
%Returns acquisition time as a string

%Likely this could rewritten to a single line code...

%Sebastian Bidhult

global SET NO
if nargin<1
 no=NO;
end

sectotal = SET(NO).AcquisitionTime;

hour = sectotal/3600;
minutes = 60*(hour-floor(hour));
sec = 60*(minutes-floor(minutes));

hour = floor(hour);
minutes = floor(minutes);
sec = floor(sec);
 if length(num2str(hour))==1
     hstr=[num2str(0),num2str(hour)];
 else
     hstr=num2str(hour);
 end
 
  if length(num2str(minutes))==1
     minstr=[num2str(0),num2str(minutes)];
 else
     minstr=num2str(minutes);
 end

 if length(num2str(sec))==1
     secstr=[num2str(0),num2str(sec)];
 else
     secstr=num2str(sec);
 end

clockstr = [hstr ':' minstr ':' secstr];