%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11.  dPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
addpath('../../../../../Empirical_Data_Analysis_Code/dPCA/matlab/')
setDir;
load ([TempDatDir 'DataList.mat']);
%Define parameter grouping
% % % parameter groupings
% % % 1 - stimulus
% % % 2 - decision
% % % 3 - time
% % % [1 3] - stimulus/time interaction
% % % [2 3] - decision/time interaction
% % % [1 2] - stimulus/decision interaction
% % % [1 2 3] - rest
% % % Here we group stimulus with stimulus/time interaction etc. Don't change
% % % that if you don't know what you are doing
% % % combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% % % margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
% % % margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
margColours    = [23 100 171; 187 20 25; 150 150 150]/256;
numTrials      = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11.1.1  PCA for collected data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for nData              = [1 2 4 5]% 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat']);
%     time               = DataSetList(nData).params.timeSeries;
%     timeEvents         = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];
%     firingRates        = generateDPCAData(nDataSet, numTrials);
%     firingRatesAverage = nanmean(firingRates, ndims(firingRates));
%     X = firingRatesAverage(:,:);
%     X = bsxfun(@minus, X, mean(X,2));
%     [W,~,~] = svd(X);
%     explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
%                                     'combinedParams', combinedParams);
%     dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
%                             'time', time,...
%                             'timeEvents', timeEvents,...
%                             'explainedVar', explVar, ...
%                             'marginalizationNames', margNames, ...
%                             'marginalizationColours', margColours);
%     setPrint(40, 24, [PlotDir 'dPCA_PCA__' DataSetList(nData).name], 'pdf')
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 11.1.2  PCA in each marginalization separately for collected data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for nData            = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat']);
%     time               = DataSetList(nData).params.timeSeries;
%     timeEvents         = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];
%     firingRates        = generateDPCAData(nDataSet, numTrials);
%     firingRatesAverage = nanmean(firingRates, ndims(firingRates));
%     dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
%                         'combinedParams', combinedParams);
%     setPrint(m*6, m*4.5, [PlotDir 'dPCA_PCA_Marginalization_' DataSetList(nData).name], 'pdf')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11.1.3  dPCA without regularization for collected data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for nData              = [1 2 4 5]% 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat']);
%     time               = DataSetList(nData).params.timeSeries;
%     timeEvents         = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];
%     firingRates        = generateDPCAData(nDataSet, numTrials);
%     firingRatesAverage = nanmean(firingRates, ndims(firingRates));
%     [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
%                     'combinedParams', combinedParams);
%     explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
%                     'combinedParams', combinedParams);
%     dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%                     'explainedVar', explVar, ...
%                     'marginalizationNames', margNames, ...
%                     'marginalizationColours', margColours, ...
%                     'whichMarg', whichMarg,                 ...
%                     'time', time,                        ...
%                     'timeEvents', timeEvents,               ...
%                     'timeMarginalization', 3, ...
%                     'legendSubplot', 16);
%     setPrint(40, 24, [PlotDir 'dPCA_dPCA_NR_' DataSetList(nData).name], 'pdf')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11.1.4  dPCA with optimal regularization for collected data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numComps             = 5;

for nData            = [1 3 4 5]% 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    time               = DataSetList(nData).params.timeSeries;
    timeEvents         = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    trialNum           = ones(size(firingRatesAverage, 1), size(firingRatesAverage, 2))*numTrials;
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, ...
                trialNum, ...
                'numComps', 5, ....
                'combinedParams', combinedParams, ...
                'numRep', 10, ...  % increase this number to ~10 for better accuracy
                'filename', 'tmp_optimalLambdas.mat',...
                'display','no');
%     setPrint(12, 9, [PlotDir 'dPCA_dPCA_ORCV_' DataSetList(nData).name], 'pdf')
    [W,V,whichMarg] = dpca(firingRatesAverage, numComps, ...
                'combinedParams', combinedParams, ...
                'lambda', optimalLambda);

    explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
                'combinedParams', combinedParams);

    dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
                'explainedVar', explVar, ...
                'marginalizationNames', margNames, ...
                'marginalizationColours', margColours, ...
                'whichMarg', whichMarg,                 ...
                'time', time,                        ...
                'timeEvents', timeEvents,               ...
                'timeMarginalization', 3,           ...
                'legendSubplot', 16);

    setPrint(40, 24, [PlotDir 'dPCA_dPCA_OR_' DataSetList(nData).name], 'pdf')
    
%     S = 2;
%     % combinedParams = {{1, [1 2]}, {2}};
%     decodingClasses = {1:S, []}; % only one decision
% 
%     accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
%         'lambda', optimalLambda, ... % 
%         'numComps', 6, ...
%         'combinedParams', combinedParams, ...
%         'decodingClasses', decodingClasses, ...
%         'numRep', 100, ...        % increase to 100
%         'filename', 'tmp_classification_accuracy.mat');
% 
% 
%     accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
%         'lambda', optimalLambda, ...
%         'combinedParams', combinedParams, ...
%         'decodingClasses', decodingClasses, ...
%         'numRep', 5, ...        % increase to 100
%         'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
%         'filename', 'tmp_classification_accuracy.mat');
% 
%     componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
% 
%     dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%         'explainedVar', explVar, ...
%         'marginalizationNames', margNames, ...
%         'marginalizationColours', margColours, ...
%         'whichMarg', whichMarg,                 ...
%         'time', time,                        ...
%         'timeEvents', timeEvents,               ...
%         'timeMarginalization', 3,           ...
%         'legendSubplot', 16,                ...
%         'componentsSignif', componentsSignif);
    
end

%% Decoding


