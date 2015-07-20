
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1
%
addpath('../Func');
setDir;
ylabels = {'Firing rate(Hz)', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};

load ([TempDatDir 'DataList.mat']);
for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityExampleTrace(nDataSet, 1:length(nDataSet), DataSetList(nData).params, ylabels{nData}); 
    m                     = ceil(sqrt(length(nDataSet)));
    setPrint(m*6, m*4.5, [PlotDir 'Single_Units__' DataSetList(nData).name], 'pdf')
end