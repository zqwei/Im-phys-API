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
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Comparison_V1_1

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
    
    % spike counts
    spkCountTime         = psthTime > wholeTrialStartingTime & psthTime < wholeTrialEndTime;
    spkCountDur          = wholeTrialEndTime - wholeTrialStartingTime;
    spk_count_yes        = cell2mat(arrayfun(@(x) getSpkCount(x, spkCountTime, spkCountDur, 'unit_yes_trial'),...
                                            nDataSet, 'UniformOutput', false))*params.binsize;
    spk_count_no         = cell2mat(arrayfun(@(x) getSpkCount(x, spkCountTime, spkCountDur, 'unit_no_trial'),...
                                            nDataSet, 'UniformOutput', false))*params.binsize;    
    
    % selectivity index for period
    periodIndex(1)       = sum(psthTime<wholeTrialStartingTime);
    periodIndex(2)       = sum(psthTime<params.polein);
    periodIndex(3)       = sum(psthTime<params.poleout);
    periodIndex(4)       = sum(psthTime<0);
    periodIndex(5)       = sum(psthTime<wholeTrialEndTime);
    sigSelective         = cell2mat(arrayfun(@(x) getSelective(x, periodIndex),...
                                            nDataSet, 'UniformOutput', false));    
                                        
                                        
    preferenceIndex      = sign(spk_count_yes - spk_count_no) == 1; 
    % preferenceIndex    == 1 contra
    % preferenceIndex    == 0 ipsi
    
    cellType             = zeros(length(preferenceIndex),1);
    
    cellType( (sigSelective(:,1)|sigSelective(:,2))  & (~sigSelective(:,3)) ) = 1;
    cellType( (sigSelective(:,1)|sigSelective(:,2))  &  sigSelective(:,3))    = 2;
    cellType( ~(sigSelective(:,1)|sigSelective(:,2))  & sigSelective(:,3))    = 3;
    
    FR_diff              = bsxfun(@rdivide, (PSTH_yes_cue_aligned - PSTH_no_cue_aligned),max(abs(PSTH_yes_cue_aligned - PSTH_no_cue_aligned), [], 2));
    
    contraCellIndex      = [];
    contraCellIndex      = [contraCellIndex; find(preferenceIndex & cellType == 1)];
    contraCellIndex      = [contraCellIndex; find(preferenceIndex & cellType == 2)];
    contraCellIndex      = [contraCellIndex; find(preferenceIndex & cellType == 3)];
    
    ipsiCellIndex        = [];
    ipsiCellIndex        = [ipsiCellIndex; find(~preferenceIndex & cellType == 1)];
    ipsiCellIndex        = [ipsiCellIndex; find(~preferenceIndex & cellType == 2)];
    ipsiCellIndex        = [ipsiCellIndex; find(~preferenceIndex & cellType == 3)];
    
    figure;
    subplot(2, 10, 1:9)
    hold on;
    imagesc(psthTime, 1:length(contraCellIndex), FR_diff(contraCellIndex,:));
    plot([params.polein params.polein],[length(contraCellIndex) 1],'--k')
    plot([params.poleout params.poleout],[length(contraCellIndex) 1],'--k')
    plot([0 0],[1 length(contraCellIndex)],'--k')
    hold off
    caxis([-1 1])
    colormap(french(128, 2));
    axis([psthTime(1) psthTime(end) 1  length(contraCellIndex)])
    xlabel('Time (s)')
    ylabel('Neuron index')
    axis ij
    colorbar
    box off
    freezeColors
    subplot(2, 10, 10)
    imagesc(psthTime, 1:length(contraCellIndex), cellType(contraCellIndex,:));
    caxis([1 3])
    colormap(french(128, 2));
    axis off
    axis ij
    axis([0.5 1.5 1  length(contraCellIndex)])
    freezeColors
    
    subplot(2, 10, 11:19)
    hold on;
    imagesc(psthTime, 1:length(ipsiCellIndex), FR_diff(ipsiCellIndex,:));
    plot([params.polein params.polein],[length(ipsiCellIndex) 1],'--k')
    plot([params.poleout params.poleout],[length(ipsiCellIndex) 1],'--k')
    plot([0 0],[1 length(ipsiCellIndex)],'--k')
    axis ij
    hold off
    caxis([-1 1])
    colormap(french);
    colorbar
    axis([psthTime(1) psthTime(end) 1  length(ipsiCellIndex)])
    xlabel('Time (s)')
    ylabel('Neuron index')
    box off
    freezeColors
    subplot(2, 10, 20)
    imagesc(psthTime, 1:length(ipsiCellIndex), cellType(ipsiCellIndex,:));
    caxis([1 3])
    colormap(french);
    axis ij
    axis off
    axis([0.5 1.5 1  length(ipsiCellIndex)])
    freezeColors
    
    setPrint(16, 12, 'reproduce_Li_Fig_2c','tif')
    
    tempSelective        = cell2mat(arrayfun(@(x) getTempSelective(x),...
                                            nDataSet, 'UniformOutput', false));    

    
    tempSelective        = tempSelective.*sign(FR_diff);                                    
    % -log_p_value plot
    
    figure;
    hold on;
    imagesc(psthTime, 1:length([contraCellIndex; ipsiCellIndex]), tempSelective([contraCellIndex; ipsiCellIndex],:));
    plot([params.polein params.polein],[length([contraCellIndex; ipsiCellIndex]) 1],'--k')
    plot([params.poleout params.poleout],[length([contraCellIndex; ipsiCellIndex]) 1],'--k')
    plot([0 0],[1 length([contraCellIndex; ipsiCellIndex])],'--k')
    axis ij 
    caxis([-5 5])
    colormap(french)
    axis([psthTime(1) psthTime(end) 1  length([contraCellIndex; ipsiCellIndex])])
    setPrint(16, 12, 'minus_log_p_value_as_a_function_time','tif')
end


function psth         = getPSTH(nData, trialType)
    averagedBinData   = mean(nData.(trialType), 1);
    boxCarWindowLength= 200; % ms
    boxCarWindow      = ones(1,boxCarWindowLength)/(boxCarWindowLength/1000);
    psth              = conv(averagedBinData, boxCarWindow, 'same');
    psth              = psth(201:end-200);
end

function spkCount     = getSpkCount(nData, spkCountTime, spkCountDur, trialType)
    averagedBinData   = mean(nData.(trialType), 1);
    spkCount          = sum(averagedBinData(spkCountTime))/spkCountDur;
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