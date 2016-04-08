% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0
%
% Comparison list
%
% 1. Reproduce Figures in Nuo's paper (Li et al., 2015, Nature)
% 2. PSTH time from Nuo: -3.3 to 1.8
%    Epoch time from Nuo: -3.1, -2.6, -1.4, 0, 1.3
% 3. My sorting algorithm
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Comparison_V1_1_1

    addpath('../Func');
    setDir;
    load ([TempDatDir 'LiAnalysis_DataList.mat']);
    
    for nData           = 1:length(DataSetList)
        load([TempDatDir DataSetList(nData).name '.mat']);
        cellType        = [DataSetList(nData).cellinfo(:).cellType];
        pryCellIndex    = cellType == 1;
        nDataSet        = nDataSet(pryCellIndex);
        PSTHStartTime   = -3.5;
        PSTHEndTime     = 2.0;
        wholeTrialStartingTime = -3.1; % pole in = -2.6 % pole out = -1.4
        wholeTrialEndTime      = 1.3; % these two values are from Nuo
        plotSelectivityIndex(nDataSet, PSTHStartTime, PSTHEndTime, wholeTrialStartingTime, wholeTrialEndTime, DataSetList(nData).params);
    end
    
end


function plotSelectivityIndex (nDataSet, PSTHStartTime, PSTHEndTime, wholeTrialStartingTime, wholeTrialEndTime, params)
    
    % psth
    psthTime             = PSTHStartTime:.001:PSTHEndTime;
    psthTime             = psthTime(201:end-200);
    PSTH_yes_cue_aligned = cell2mat(arrayfun(@(x) getPSTH(x, 'unit_yes_trial'),...
                                            nDataSet, 'UniformOutput', false))*params.binsize;
    PSTH_no_cue_aligned  = cell2mat(arrayfun(@(x) getPSTH(x, 'unit_no_trial'),...
                                            nDataSet, 'UniformOutput', false))*params.binsize;
                                           
    FR_diff              = bsxfun(@rdivide, (PSTH_yes_cue_aligned - PSTH_no_cue_aligned),max(abs(PSTH_yes_cue_aligned - PSTH_no_cue_aligned), [], 2));
    
    
    tempSelective        = cell2mat(arrayfun(@(x) getTempSelective(x),...
                                            nDataSet, 'UniformOutput', false));  
                                        
    periodIndex(1)       = sum(psthTime<wholeTrialStartingTime);
    periodIndex(2)       = sum(psthTime<params.polein);
    periodIndex(3)       = sum(psthTime<params.poleout);
    periodIndex(4)       = sum(psthTime<0);
    periodIndex(5)       = sum(psthTime<wholeTrialEndTime);
    sigSelective         = cell2mat(arrayfun(@(x) getSelective(x, periodIndex),...
                                            nDataSet, 'UniformOutput', false));                                        
    selctiveNeurons      = sum(sigSelective, 2) > 0;
    
    tempSelective        = tempSelective(selctiveNeurons, :);
    FR_diff              = FR_diff(selctiveNeurons, :);
    
    actMat                    = tempSelective;
    zScores                   = tempSelective.*sign(FR_diff);                                    
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
    similaritySort                = similaritySort(bumpStartPoint(similaritySort)<numT);
    
    figure;
    hold on;
    imagesc(psthTime, 1:length(similaritySort), zScores(similaritySort,:));
    plot([params.polein params.polein],[length(similaritySort) 1],'--k')
    plot([params.poleout params.poleout],[length(similaritySort) 1],'--k')
    plot([0 0],[1 length(similaritySort)],'--k')
    axis ij 
    caxis([-5 5])
    colormap(french)
    axis([psthTime(1) psthTime(end) 1  length(similaritySort)])
    setPrint(16, 12, 'minus_log_p_value_as_a_function_time_new_soring','tif')
end


function psth         = getPSTH(nData, trialType)
    averagedBinData   = mean(nData.(trialType), 1);
    boxCarWindowLength= 200; % ms
    boxCarWindow      = ones(1,boxCarWindowLength)/(boxCarWindowLength/1000);
    psth              = conv(averagedBinData, boxCarWindow, 'same');
    psth              = psth(201:end-200);
end


function tempSelective= getTempSelective(nData)
    boxCarWindowLength= 200; % ms
    boxCarWindow      = ones(1,boxCarWindowLength)/(boxCarWindowLength/1000);
    yesActMat         = nData.unit_yes_trial;
    noActMat          = nData.unit_no_trial;
    
    for nAct          = 1:size(yesActMat, 1)
        yesActMat(nAct,:) = conv(yesActMat(nAct,:), boxCarWindow, 'same');
    end
    
    for nAct          = 1:size(noActMat, 1)
        noActMat(nAct,:)  = conv(noActMat(nAct,:), boxCarWindow, 'same');
    end
    
    yesActMat         = yesActMat(:,201:end-200);
    noActMat          = noActMat(:,201:end-200);
    
    [~, tempSelective, ~] = ttest2(yesActMat, noActMat,'dim',1);
    tempSelective(~isnan(tempSelective)) = -log(tempSelective(~isnan(tempSelective)));
    tempSelective(isnan(tempSelective))  = 0;

end


function sigSelective = getSelective(nData, periodIndex)
    sigSelective      = false(1,length(periodIndex) - 2);
    for nPeriod       = 2:length(periodIndex)-1
        yesData       = mean(nData.unit_yes_trial(:,periodIndex(nPeriod):periodIndex(nPeriod+1)),2);
        noData        = mean(nData.unit_no_trial(:,periodIndex(nPeriod):periodIndex(nPeriod+1)),2);
        h             = ttest2(yesData, noData);
        if isnan(h); h = false; end
        sigSelective(nPeriod-1) = h;
    end
end