function plotCaSelectivityIndex (nDataSet, wholeTrialStartingTime, wholeTrialEndTime, params)
    
    % psth
    psthTime             = params.timeSeries;
    PSTH_yes_cue_aligned = cell2mat(arrayfun(@(x) mean(x.unit_yes_trial),...
                                            nDataSet, 'UniformOutput', false))*params.binsize;
    PSTH_no_cue_aligned  = cell2mat(arrayfun(@(x) mean(x.unit_no_trial),...
                                            nDataSet, 'UniformOutput', false))*params.binsize;
    
    % spike counts
    % spkCountTime         = psthTime > wholeTrialStartingTime & psthTime < wholeTrialEndTime;
    spkCountTime         = psthTime > params.polein & psthTime < wholeTrialEndTime;
    % spkCountDur          = wholeTrialEndTime - wholeTrialStartingTime;
    spkCountDur          = wholeTrialEndTime - params.polein;
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
                                        
    % based on the spike counts across all three trial epochs                                    
    preferenceIndex      = sign(spk_count_yes - spk_count_no); 
    % preferenceIndex    == 1 contra
    % preferenceIndex    == 0 ipsi
    
    cellType             = zeros(length(preferenceIndex),1);
    
    cellType( (sigSelective(:,1)|sigSelective(:,2))  & (~sigSelective(:,3)) ) = 1;
    cellType( (sigSelective(:,1)|sigSelective(:,2))  &  sigSelective(:,3))    = 2;
    cellType( ~(sigSelective(:,1)|sigSelective(:,2))  & sigSelective(:,3))    = 3;
    
    FR_diff              = bsxfun(@rdivide, (PSTH_yes_cue_aligned - PSTH_no_cue_aligned),max(abs(PSTH_yes_cue_aligned - PSTH_no_cue_aligned), [], 2));
    
    contraCellIndex      = [];
    contraCellIndex      = [contraCellIndex; find(preferenceIndex==1 & cellType == 1)];
    contraCellIndex      = [contraCellIndex; find(preferenceIndex==1 & cellType == 2)];
    contraCellIndex      = [contraCellIndex; find(preferenceIndex==1 & cellType == 3)];
    
    ipsiCellIndex        = [];
    ipsiCellIndex        = [ipsiCellIndex; find(preferenceIndex==-1 & cellType == 1)];
    ipsiCellIndex        = [ipsiCellIndex; find(preferenceIndex==-1 & cellType == 2)];
    ipsiCellIndex        = [ipsiCellIndex; find(preferenceIndex==-1 & cellType == 3)];
    
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
    box off
    set(gca, 'YTick', [1  length(contraCellIndex)])
    set(gca, 'TickDir', 'out') 
    freezeColors
    
    
    subplot(2, 10, 10)
    imagesc(psthTime, 1:length(contraCellIndex), cellType(contraCellIndex,:));
    caxis([1 4])
%     colormap(gray);
    colormap(parula);
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
    axis([psthTime(1) psthTime(end) 1  length(ipsiCellIndex)])
    xlabel('Time (s)')
    ylabel('Neuron index')
    set(gca, 'YTick', [1  length(ipsiCellIndex)])
    set(gca, 'TickDir', 'out') 
    box off
    freezeColors
    
    subplot(2, 10, 20)
    imagesc(psthTime, 1:length(ipsiCellIndex), cellType(ipsiCellIndex,:));
    caxis([1 4])
%     colormap(gray);
    colormap(parula);
    axis ij
    axis off
    axis([0.5 1.5 1  length(ipsiCellIndex)])    
    freezeColors

end