%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first p-value ramping time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);
% Gaussian filter for spiking data
sigma                         = 0.1 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);


for nData                     = [3 4]
    load([TempDatDir DataSetList(nData).name '.mat']);
    numUnits                  = length(nDataSet);
    pValue                    = applyFuncToCompareTrialType(nDataSet, @pValueTTest2);
    meanDiffValue             = applyFuncToCompareTrialType(nDataSet, @meanDiff);
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
    [~, similaritySort]           = sortrows([bumpStartPoint, bumpSize, bumpSign], [-3 -1 -2]);    

    h = figure;
    hold on
    imagesc(DataSetList(nData).params.timeSeries, 1:numUnits, zScores(similaritySort,:));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)])
    ylim([1 numUnits])
    set(gca, 'YTick', [1 numUnits])
    colormap(h, french(128,2))
%     colorbar;
    caxis([-5 5])
    axis xy
    xlabel('Time (s)')
    ylabel('Neuron Index')
    box off;
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsFPVT/SingleUnitsZScore_' DataSetList(nData).name])
end


close all