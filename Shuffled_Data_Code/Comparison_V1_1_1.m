%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Raw activity of all neurons -- neuron classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1_1
%
%
%
% Neurons are classified into 4 groups
%
% group 1. stimulus driven neurons
% group 2. response driven neurons
% group 3. preparatory neurons
% group 4. non-selective neurons



addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

skipPlot                = false;
if ~skipPlot
    for nData             = 1:length(DataSetList)
        load([TempDatDir DataSetList(nData).name '.mat'])
        plotMeanActivityImagescWithCellClass(nDataSet, DataSetList(nData).params, DataSetList(nData).cellinfo, [], []);
%        m                     = ceil(sqrt(length(nDataSet)));
%        setPrint(m*4, m*3, ['Plot/Single_Units_' DataSetList(nData).name], 'pdf')
    end
end


