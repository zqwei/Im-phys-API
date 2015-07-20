%
% Get the data from the data files of Nuo and Tsai-Wen
%
% 


addpath('../Func');
setDir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDataSet               = compileSpikeData(SpikingDataDir, SpikeFileList); %#ok<NASGU>
DataSetList(1).name    = 'Compile_Spikes';
save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow virus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDataSet               = compileCaImagingData(CaImagingShortDelaySlowVirusDir, CaImagingShortDelaySlowVirusFileList);
DataSetList(2).name    = 'Compile_Ca_Slow_Short_Delay_Virus';
save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet');