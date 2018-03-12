addpath('../Func');
setDir;

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
load([TempDatDir DataSetList(1).name '.mat'])
params    = DataSetList(1).params;
binSize   = params.binsize;

numTime   = length(params.timeSeries);
numNeuron = length(nDataSet);
timePoints= timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

% instanteous noise level
var_mat   = nan(numNeuron, numTime * 2); 
sig_mat   = nan(numNeuron, numTime * 2); 

for nNeuron = 1:numNeuron
   sig_mat(nNeuron, 1:numTime)     = mean(nDataSet(nNeuron).unit_yes_trial);
   var_mat(nNeuron, 1:numTime)     = var(nDataSet(nNeuron).unit_yes_trial);
   sig_mat(nNeuron, 1+numTime:end) = mean(nDataSet(nNeuron).unit_no_trial);
   var_mat(nNeuron, 1+numTime:end) = var(nDataSet(nNeuron).unit_no_trial);
end

noise_dist = (var_mat(:)*params.binsize^2)./(sig_mat(:)*params.binsize);
per_list   = 2:98;
noise_factor_list = prctile(noise_dist, per_list);

params            = DataSetList(2).params;
indexDatasets = 4;% [3, 4];
% 2: short Ca GP517
% 3: short Ca slow
% 4: short Ca slow virus
% 5: long Ca fast
% 6: long Ca slow

timeTag        = 8:60;
numComps       = 3;
combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 1;
numTrialsLDA   = 500;
numTestTrials  = 200;
numTrainingTrials = numTrialsLDA - numTestTrials;
numRandPickUnits  = 50;


factorSet         = [0 0 2.5 5.5 0 0];

for nData         = indexDatasets    
    
    load([TempDatDir 'directDeconv_' DataSetList(nData).name '.mat'], 'yesData', 'noData');   
    analysisMat   = repmat(struct('nParaSet',1, 'pCa0', 1, ...
                                'pTaud', 1, 'sizeGroup', 1),900, 1);
    numUnit       = size(yesData, 2);
    allYesData    = yesData;
    allNoData     = noData;
    
    randTau       = randi(97, 1000, numUnit);
    randFF        = randi(97, 1000, numUnit);
    nDataSet      = repmat(struct('unit_yes_trial', 1, 'unit_no_trial', 1),numUnit, 1);
    validUnit     = true(numUnit, 1);
    
    for nParaSet  = 1:1000        
        
        for nUnit = 1:numUnit
            ffactor = noise_factor_list(randFF(nParaSet, nUnit));
            mean_yes_trial = squeeze(allYesData(randTau(nParaSet, nUnit), nUnit, :))';
            std_yes_trial  = ones(numTrialsLDA/2, 1) * sqrt(mean_yes_trial* ffactor/params.binsize) .* randn(numTrialsLDA/2, 77);
            min_std_yes_trial = min(std_yes_trial)./mean_yes_trial;
            mean_no_trial  = squeeze(allNoData(randTau(nParaSet, nUnit), nUnit, :))';
            std_no_trial   = ones(numTrialsLDA/2, 1) * sqrt(mean_no_trial* ffactor/params.binsize) .* randn(numTrialsLDA/2, 77);
            min_std_no_trial = min(std_no_trial)./mean_no_trial;
            nfactor          = factorSet(nData);
            nDataSet(nUnit).unit_yes_trial   = mean_yes_trial * nfactor + std_yes_trial*sqrt(nfactor);
            nDataSet(nUnit).unit_no_trial    = mean_no_trial * nfactor + std_no_trial*sqrt(nfactor);
            if sum(isnan(nDataSet(nUnit).unit_yes_trial(:)))+sum(isnan(nDataSet(nUnit).unit_no_trial(:)))>0
                validUnit(nUnit) = false;
            end
        end
        
        nDataSet  = nDataSet(validUnit);
        if size(nDataSet, 1) == 1
            nDataSet = nDataSet';
        end
        
        % cell type data
        unitGroup = getLogPValueTscoreC2STime(nDataSet, params);       
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
            decodability(nFold,:) = decodabilityLDA(nSessionData, trainingTargets, testTargets);
        end
        analysisMat(nParaSet).nParaSet     = nParaSet;
        analysisMat(nParaSet).sizeGroup    = sizeGroup;
        analysisMat(nParaSet).peakiness    = peakiness;
        analysisMat(nParaSet).PCAVar       = PCAVar;
        analysisMat(nParaSet).decodability = decodability;
    end
    
    save([TempDatDir 'ResultsCompiledC2S_' DataSetList(nData).name '.mat'], 'analysisMat', 'randTau', 'randFF');
end