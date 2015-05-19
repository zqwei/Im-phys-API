%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Imagec for the different dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotMeanActivityImagesc(SpikeDataSet, paramsSpike.timeSeries)
% plotMeanActivityImagesc(FakeCaImagingDataSet, paramsSpike.timeSeries)
% plotMeanActivityImagesc(CaImagingShortDelay, paramsROIS.timeWindowIndexRange / frameRate)
% plotMeanActivityImagesc(CaImagingLongDelay, paramsROIL.timeWindowIndexRange / frameRate)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 1.1
% for loop for plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load ('TempDat/DataList.mat');
% for nData             = 1:length(DataSetList)
%     load(['TempDat/' DataSetList(nData).name '.mat'])
%     plotMeanActivityImagesc(nDataSet, DataSetList(nData).params);  
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Neural Traces for example neurons (randomly chosen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampleSize            = 16;
% m                     = ceil(sqrt(sampleSize));
% numUnits              = 200;
% sampleSeq             = randperm(numUnits);
% sampleSeq             = sampleSeq(1:sampleSize);
% % plotMeanActivityExampleTrace(SpikeDataSet, paramsSpike.timeSeries, sampleSize, -2.6, -1.3)
% % plotMeanActivityExampleTrace(FakeCaImagingDataSet, paramsSpike.timeSeries, sampleSize, -2.6, -1.3)
% % plotMeanActivityExampleTrace(CaImagingShortDelay, paramsROIS.timeWindowIndexRange / frameRate, sampleSize, -2.6, -1.4)
% % plotMeanActivityExampleTrace(CaImagingLongDelay, paramsROIL.timeWindowIndexRange / frameRate, sampleSize, -4.2, -3.0)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 1.1
% for loop for plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampleSize            = 16;
% m                     = ceil(sqrt(sampleSize));
% numUnits              = 200;
% sampleSeq             = randperm(numUnits);
% sampleSeq             = sampleSeq(1:sampleSize);
% load ('TempDat/DataList.mat');
% for nData             = 1:length(DataSetList)
%     load(['TempDat/' DataSetList(nData).name '.mat'])
%     plotMeanActivityExampleTrace(nDataSet, sampleSeq, DataSetList(nData).params); 
%     setPrint(m*4, m*3, ['Plot/Single_Units_' DataSetList(nData).name], 'png')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Neural Traces for all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skipPlot                = true;
if ~skipPlot
    for nData             = fileToAnalysis
        load([TempDatDir DataSetList(nData).name '.mat'])
        plotMeanActivityExampleTrace(nDataSet, 1:length(nDataSet), DataSetList(nData).params); 
        m                     = ceil(sqrt(length(nDataSet)));
        setPrint(m*4, m*3, ['Plot/Single_Units_' DataSetList(nData).name], 'pdf')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 2.0
% based on Shaul's code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};
% yAxes_set                = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0];
% lowFiringThres          = [20, 0.5, 0.5, 0.5, 0.5];
% skipPlot                = false;
% if ~skipPlot
%     for nData             = 1:1%2:length(DataSetList)
%         load([TempDatDir DataSetList(nData).name '.mat'])
%         plotMeanActivityImagescWithSort(nDataSet, DataSetList(nData).params, [], 0, ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%         setPrint(6*4, 3*3, [PlotDir '/Single_UnitsImagescWithSort_' DataSetList(nData).name], 'tif')
%     end
% end



ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F' , 'dF/F'};
yAxes_set               = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0 ; -0.5 2.0];
lowFiringThres          = [15, 0.3, 0.3, 0.3, 0.3, 0.3];
skipPlot                = true;
if ~skipPlot
    if ~exist([PlotDir '/Single_UnitsImagescWithSort'],'dir')
        mkdir([PlotDir '/Single_UnitsImagescWithSort'])
    end
    for nData             = fileToAnalysis
        load([TempDatDir DataSetList(nData).name '.mat'])
        plotMeanActivityImagescWithSortWithCellinfo(nDataSet(ephysCellIndex{nData}), DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%         title(DataSetList(nData).name,'interpreter','none')
        setPrint(6*4, 3*3, [PlotDir '/Single_UnitsImagescWithSort/Single_UnitsImagescWithSort_' DataSetList(nData).name], 'tif')
    end
end

close all;


ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F' , 'dF/F'};
yAxes_set               = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0 ; -0.5 2.0];
lowFiringThres          = [15, 0.3, 0.3, 0.3, 0.3, 0.3];
skipPlot                = false;
if ~skipPlot
    if ~exist([PlotDir '/Single_UnitsImagescWithSort'],'dir')
        mkdir([PlotDir '/Single_UnitsImagescWithSort'])
    end
    for nData             = 1:length(DataSetList)
        load([TempDatDir DataSetList(nData).name '.mat'])
        plotMeanActivityImagescWithSortWithCellinfo(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%         title(DataSetList(nData).name,'interpreter','none')
        setPrint(6*4, 3*3, [PlotDir '/Single_UnitsImagescWithSort/Single_UnitsImagescWithSort_' DataSetList(nData).name], 'tif')
    end
end

close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 3.0
% based on 3D location information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};
% yAxes_set                = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0];
% lowFiringThres          = [20, 0.5, 0.5, 0.5, 0.5];
% skipPlot                = false;
% if ~skipPlot
%     for nData             = 1:length(DataSetList)
%         load([TempDatDir DataSetList(nData).name '.mat'])
%         plotMeanActivityImagescWithSort(nDataSet, DataSetList(nData).params, [], 0, ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%         setPrint(6*4, 3*3, [PlotDir '/Single_UnitsImagescWithSort_' DataSetList(nData).name], 'tif')
%     end
% end
