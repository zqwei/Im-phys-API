addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params            = DataSetList(1).params;

dataSetNames   = {'Modeled_6s_AAV_nRSParaSet_', 'Modeled_GP43_nRSParaSet_'};

timeTag        = 8:60;
numComps       = 3;
combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 100;
numTrialsLDA   = 500;
numTestTrials  = 200;
numTrainingTrials   = numTrialsLDA - numTestTrials;
numRandPickUnits = 50;

for nData         = 1:2
    for nParaSet  = 1:1000
        load([TempDatDir dataSetNames{nData} num2str(nParaSet, '%04d') '.mat'])    
        
        % cell type data
        unitGroup = getLogPValueTscoreSpikeTime(nDataSet, params);       
        sizeGroup = histcounts(unitGroup, 0:3); % record this for cell type data
        
        % peakiness
        numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
        yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
        noProfileMatrix     = yesProfileMatrix;
        positivePeak        = false(length(nDataSet));
        for nUnit        = 1:length(nDataSet)
            yesData      = mean(nDataSet(nUnit).unit_yes_trial);
            noData       = mean(nDataSet(nUnit).unit_no_trial);
            maxData      = max([yesData, noData]);
            minData      = min([yesData, noData]);
            rData        = (maxData - minData);
            yesData      = (yesData - minData)/(maxData - minData);
            noData       = (noData - minData)/(maxData - minData);
            yesProfileMatrix(nUnit, :)    = yesData;
            noProfileMatrix(nUnit, :)     = noData;
            positivePeak(nUnit)       = mean(yesData(1:8)) <= mean(yesData(9:47)) ...
                                       || mean(noData(1:8)) <= mean(noData(9:47));
        end
        actMat        = [yesProfileMatrix, noProfileMatrix];
        actMat        = actMat(positivePeak, :);
        [~, maxId]    = max(actMat, [], 2);
        countMaxId    = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
        peakiness     = std(countMaxId([timeTag, timeTag+77]));
        
        % pca
        evMat              = zeros(numFold, length(combinedParams), numComps);
        firingRates        = generateDPCAData(nDataSet, numTrials);
        firingRatesAverage = nanmean(firingRates, ndims(firingRates));
        pcaX               = firingRatesAverage(:,:);
        firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
        pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
        Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
        totalVar           = sum(sum(pcaX.^2));
        [~, ~, Wpca] = svd(pcaX');
        PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
        for i=1:length(Xmargs)
            PCAmargVar(i,:)= sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
        end
        PCAVar             = PCAmargVar(:, 1:numComps)';
        
        %lda
        decodability            = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));        
        for nFold    = 1:numFold
            trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
            trainingTargets     = trainingTargets(randperm(numTrainingTrials));
            testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
            testTargets         = testTargets(randperm(numTestTrials));
            totTargets          = [testTargets; trainingTargets];
            trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
            testDecisions       = testTargets(randperm(numTestTrials));
            totDecisions        = [testDecisions; trainingDecisions];
            randPickUnits       = randperm(length(nDataSet));
            randPickUnits       = randPickUnits(1:numRandPickUnits);
            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrialsLDA), trainingTargets, testTargets);
        end
        clear nDataSet
        
        save([TempDatDir 'Results_' dataSetNames{nData} num2str(nParaSet, '%04d') '.mat'], 'sizeGroup', 'countMaxId', 'peakiness', 'PCAVar', 'decodability')
    end
end