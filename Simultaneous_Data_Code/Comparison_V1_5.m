%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  Per session decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.1  Per session decodability over time -- shuffled trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_5

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);

addNoise         = [1 0 0 0 0];

numFold          = 30;



for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    numSession       = length(nDataSet3D);
    m                = ceil(numSession/4);
    numUnit          = arrayfun(@(x) length(x.nUnit), nDataSet3D, 'UniformOutput', false);
    [numUnit, sortUnit]    = sort([numUnit{:}],'descend');
    
    nDataSet3D       = nDataSet3D(sortUnit);
    figure;
    for nPlot        = 1:numSession
        subplot(m, 4, nPlot)
        nSessionDataOrginal = [permute(nDataSet3D(nPlot).unit_yes_trial,[2 1 3]); permute(nDataSet3D(nPlot).unit_no_trial,[2 1 3])];
        nTargetsOrginal     = [true(size(nDataSet3D(nPlot).unit_yes_trial, 2),1); false(size(nDataSet3D(nPlot).unit_no_trial, 2),1)];
        n_yes        = size(nDataSet3D(nPlot).unit_yes_trial, 2);
        n_no         = size(nDataSet3D(nPlot).unit_no_trial, 2);
        minTrial     = min(n_yes, n_no);
        
        if numUnit(nPlot) < minTrial - 10
        
            decodability = zeros(numFold, size(nSessionDataOrginal, 3));

            for nFold    = 1:numFold        
                randSeqSession = [randperm(n_yes) randperm(n_no)+n_yes];
                randSeqSession = [randSeqSession(1:minTrial); randSeqSession(n_yes+1:n_yes+minTrial)];
                randSeqSession = randSeqSession(:);
                nTargets     = nTargetsOrginal(randSeqSession);
                nSessionData = nSessionDataOrginal(randSeqSession, :, :);
                nTest        = nTargets(1:minTrial);
                nTrain       = nTargets(minTrial+1:end);
                % nSessionData = [nSessionData; nSessionData]; %#ok<AGROW>
                decodability(nFold, :) = decodabilityLDA(nSessionData+randn(size(nSessionData))*1e-3/sqrt(length(nTargets))*addNoise(nData), nTrain, nTest);
            end

            hold on
            plot(DataSetList(nData).params.timeSeries, mean(decodability,1), 'k', 'linewid',1);
    %         plot(DataSetList(nData).params.timeSeries, mean(decodability,1)+std(decodability,[],1)/sqrt(numFold), 'r', 'linewid',0.5);
    %         plot(DataSetList(nData).params.timeSeries, mean(decodability,1)-std(decodability,[],1)/sqrt(numFold), 'r', 'linewid',0.5);
            xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
            ylim([0.4 1])
            set(gca,'YTick',0.4:0.2:1)
            gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
            box off;
            hold off;
            xlabel('Time (s)');
            ylabel('Decodability');
            title([num2str(length(nDataSet3D(nPlot).nUnit)) ' Neurons'])
        end
    end
    setPrint(4*4, m*3, [PlotDir 'Single_Session_Decodability_' DataSetList(nData).name], 'pdf')
end



