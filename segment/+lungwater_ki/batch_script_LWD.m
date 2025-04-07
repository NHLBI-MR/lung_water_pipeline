segment
global SET NO DATA
warning('off','all');
clc

doAnalysis=1;
doLungSeg=1; %1 if you want to run U-net, 0 if you have segmentations already
doPlot=1;
doSave=1;
doExport=1; doHeader=1; output=[];
bodymaskslices=[]; doLiverRoi=0;

pathname='C:\Users\username\Documents\Data\'; %enter the folder with your segment files


files=dir([pathname,'*.mat'])
filenameXLSX='results_LWD_exportF.xlsx';
files.name

segment('filecloseall_Callback');
%%
%loop through files
for i=1:length(files)

    % --- Load file
    DATA.Silent = true;  % Turn on "silent" mode to avoid to much update on screen when loading etc.
    filename = files(i).name(1:end-4);
    disp(sprintf('Loading %s',filename));

    % Make sure a fresh start
    SET = []; NO = 1;
    load([pathname filesep files(i).name],'-mat'); % Load
    % Assign
    SET = setstruct;
    clear setstruct;
    % Call to intialize all variables correcly after loaded data.
    DATA.Preview.PreviewFile = [pathname filesep files(i).name];

    openfile('setupstacksfrommat',NO);
    segment('renderstacksfrommat');

    %Does all graphical updates
    if not(isempty(DATA.ViewMatrix))
        NO=1; DATA.init_graphics;
    end

    %update ROIs so they are consistent across timeframes
    for no=1:length(SET)
        tf_draw=1; %I AM ASSUMING YOU MADE MANUAL CORRECTIONS IN THE FIRST TIMEFRAME, CHANGE IF YOU USED ANOTHER ONE
        ind = lungwater_ki.indexROISlice(no,'Lung'); %find ROIs names lung
        for rois=ind %loop through Lung ROI's
            if sum(sum(SET(no).Roi(rois).X-SET(no).Roi(rois).X(:,tf_draw),2))~=0 %check if ROIs are changing over time
                for tf=1:size(SET(no).Roi(rois).X,2) %loop over timeframes to update ROI contours
                    SET(no).Roi(rois).X(:,tf)=SET(no).Roi(rois).X(:,tf_draw);
                    SET(no).Roi(rois).Y(:,tf)=SET(no).Roi(rois).Y(:,tf_draw);
                end
            end
        end
    end

    if doAnalysis %run LWD analysis  
        lungwater_ki.dynamic_lungwater_pipeline([1:length(SET)], doLungSeg, doPlot, bodymaskslices, doLiverRoi, filename); %runs lung water pipeline in all SETs
    end

    %save .mat-file
    if doSave
        disp('Saving segment file.');
        DATA.Silent=true;
        filemenu('saveallas_helper', pathname, filename);
        DATA.Silent=false;
        disp('File saved');
    end

    %export results
    if doExport
        out=lungwater_ki.exportDynamicLWDSegment([1:length(SET)], doHeader);
        output=[output;out];
        doHeader=0; %only print header in the first row of all subjects
    end
end

if doExport
    % xlswrite(fullfile(pathname,filenameXLSX),output);
    array2table(output)
    writecell(output, fullfile(pathname,filenameXLSX));
    % type(filenameXLSX);
end

DATA.Silent = false; %turn off silent mode
disp('Finished');
