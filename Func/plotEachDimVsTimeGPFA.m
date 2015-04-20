%
%
%
% Ziqiang Wei
% 

function plotEachDimVsTimeGPFA(seq, params)

    
    nCols     = 4;

    f = figure;
    pos = get(gcf, 'position');
    set(f, 'position', [pos(1) pos(2) 2*pos(3) pos(4)]);

    Xall = [seq.xorth];
    xMax = ceil(10 * max(abs(Xall(:)))) / 10; % round max value to next highest 1e-1

    ytk       = [-xMax 0 xMax];
    nRows     = ceil(size(Xall, 1) / nCols);

    yes_trial  = [seq(:).trialType];
    [xDim, nT] = size(seq(1).xorth);

    xorths     = zeros(length(seq), xDim, nT);

    for  n     = 1:length(seq)
      xorths(n, :, :)  = seq(n).xorth;
    end
  
  
        
    for k = 1:xDim
        subplot(nRows, nCols, k);
        hold on;
        shadedErrorBar(params.timeSeries, mean(xorths(yes_trial==1,k,:),1), std(xorths(yes_trial==1,k,:),[],1)/sqrt(sum(yes_trial==1)), {'-b', 'linewid', 1.0}, 0.5);
        shadedErrorBar(params.timeSeries, mean(xorths(yes_trial==0,k,:),1), std(xorths(yes_trial==0,k,:),[],1)/sqrt(sum(yes_trial==0)), {'-r', 'linewid', 1.0}, 0.5);
        axis([params.timeSeries(1) params.timeSeries(end) 1.1*min(ytk) 1.1*max(ytk)]);
        gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
        str = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',k);
        title(str, 'interpreter', 'latex', 'fontsize', 16);
        xlabel('Time (ms)');
    end