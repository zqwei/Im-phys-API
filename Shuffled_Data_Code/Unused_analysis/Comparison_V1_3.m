%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Selectivity over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_3
%
% Selectivity is computed as two sample z-score in two trial type condition

% plotZScoreImagesc (SpikeDataSet, paramsSpike.timeSeries)
% plotZScoreImagesc (FakeCaImagingDataSet, paramsSpike.timeSeries)
% % For FakeCaImagingDataSet, its z-score spans in a large range.
% plotZScoreImagesc (CaImagingShortDelay, paramsROIS.timeWindowIndexRange / frameRate)
% plotZScoreImagesc (CaImagingLongDelay, paramsROIL.timeWindowIndexRange / frameRate)

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotZScoreImagesc(nDataSet, DataSetList(nData).params); 
    setPrint(8, 6, [PlotDir 'Single_Units_zScore_' DataSetList(nData).name], 'tif')
end