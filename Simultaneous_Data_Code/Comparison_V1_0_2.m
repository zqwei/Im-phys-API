% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.7
% 

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
% Gaussian filter for spiking data
sigma                         = 0.10 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 30;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

if ~exist([PlotDir '/All_Data_Summary'],'dir')
    mkdir([PlotDir '/All_Data_Summary'])
end

figure;
hold on;
legendString                  = cell(length(fileToAnalysis),1);
for nData                     = fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat']);
    numUnits                  = length(nDataSet);
    pValue                    = applyFuncToCompareTrialType(nDataSet, @pValueTTest2, filterInUse);
    meanDiffValue             = applyFuncToCompareTrialType(nDataSet, @meanDiff, filterInUse);
    logPValue                 = -log(pValue);
    % yes  -- blue trial
    % no   -- red trial
    zScores                   = -sign(meanDiffValue).*logPValue;     
    actMat                    = logPValue;
    numT                      = size(zScores,2);
    bumpActThres              = 3; % > bumpActThres considering as a bump % 3 = -log(0.05)
    bumpMat                   = actMat > bumpActThres;
    bumpSize                  = ones(size(actMat,1),1);
    bumpStartPoint            = nan(size(actMat,1),1);
    bumpSign                  = ones(size(actMat,1),1);

    for nUnit = 1: size(actMat,1)
        diffActMat                = diff([bumpMat(nUnit,:),0]);
        beginPoint                = find(diffActMat==1);
        endPoint                  = find(diffActMat==-1);
        if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
        [bumplength, bumpIndex]   = max(endPoint - beginPoint);
        if isempty(bumplength); bumplength = 0; end
        bumpSize(nUnit)           = bumplength;
        if bumplength==0
            bumpStartPoint(nUnit) = numT;
            bumpSign(nUnit)       = -1;
        else 
            bumpStartPoint(nUnit) = beginPoint(bumpIndex);
            bumpSign(nUnit)       = actMat(nUnit, bumpStartPoint(nUnit))/zScores(nUnit, bumpStartPoint(nUnit));
        end
    end

    startTimePoint               = DataSetList(nData).params.timeSeries(bumpStartPoint);    
    nCell                        = histcounts(startTimePoint, DataSetList(nData).params.timeSeries(1:3:end-3));
    plot(DataSetList(nData).params.timeSeries(4:3:end-3), cumsum(nCell)/sum(nCell));
    legendString(fileToAnalysis==nData) = {[strrep(DataSetList(nData).name,'_',' '), ': ', num2str(sum(nCell)) ' cells']}; 
end
legend(legendString, 'Location', 'northwest');
legend('boxoff')
xlim([DataSetList(min(nData,4)).params.timeSeries(1) DataSetList(min(nData,4)).params.timeSeries(end)])
ylim([0 1])
ylabel('Fraction of Cell')
xlabel('First Sig. P Val. Time (s)')
setPrint(8, 6, [PlotDir 'All_Data_Summary/Cropped_SimultaneousData_FPVT'], 'pdf')