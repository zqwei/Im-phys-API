    
% %     sigma                         = 0.05 / params.binsize; % 200 ms
% %     filterLength                  = 11;
% %     filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
% %     filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
% %     filterInUse                   = filterInUse / sum (filterInUse); 
% %     
% %     yesActMat                     = getGaussianPSTH (filterInUse, yesActMat, 2);
% %     noActMat                      = getGaussianPSTH (filterInUse, noActMat, 2);
%     
%     otherIndex = unitGroup~=0 & cellType == 1 & validDepth';
%     
%     figure;
%     subplot(1, 2, 1)
%     hold on
%     shadedErrorBar(params.timeSeries, mean(noActMat(contraIndex & otherIndex,:)), sem(noActMat(contraIndex& otherIndex,:))/2,'-r')
%     shadedErrorBar(params.timeSeries, mean(yesActMat(contraIndex& otherIndex,:)), sem(yesActMat(contraIndex& otherIndex,:))/2,'-b')
%     ylim([2, 10])
%     title('contra neuron only')
%     gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%     xlim([params.timeSeries(2) params.timeSeries(end-1)])
%     ylim([2, 10])
%     xlabel('Time (ms)')
%     box off
%     ylabel('Mean activity')
%     
%     baseline1 = mean(mean(noActMat(contraIndex, 1:8)));
%     baseline2 = mean(mean(yesActMat(contraIndex, 1:8)));
%     baseline  = (baseline1 + baseline2)/2;
%     
%     disp(max(mean(noActMat(contraIndex,:))) - baseline)
%     disp(min(mean(yesActMat(contraIndex,:))) - baseline)
% 
%     subplot(1, 2, 2)
%     hold on
%     shadedErrorBar(params.timeSeries, mean(noActMat(~contraIndex& otherIndex,:)), sem(noActMat(~contraIndex& otherIndex,:))/2,'-r')
%     shadedErrorBar(params.timeSeries, mean(yesActMat(~contraIndex& otherIndex,:)), sem(yesActMat(~contraIndex& otherIndex,:))/2,'-b')
%     ylim([2, 10])
%     title('ipsi neuron only')
%     gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%     xlim([params.timeSeries(2) params.timeSeries(end-1)])
%     ylim([2, 10])
%     xlabel('Time (ms)')
%     box off
%     ylabel('Mean activity')
%     
%     baseline1 = mean(mean(noActMat(~contraIndex, 1:8)));
%     baseline2 = mean(mean(yesActMat(~contraIndex, 1:8)));
%     baseline  = (baseline1 + baseline2)/2;
%     
%     disp(min(mean(noActMat(~contraIndex,:))) - baseline)
%     disp(max(mean(yesActMat(~contraIndex,:))) - baseline)
%     
%     setPrint(8*2, 6, [PlotDir 'SingleUnitsContraIpsi\' DataSetList(nData).name])