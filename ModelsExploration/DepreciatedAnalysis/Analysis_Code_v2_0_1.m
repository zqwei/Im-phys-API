%
% Comparison between linear model and TWC model
%
% Plot
% x -- ca++ trace
% y -- modeled trace
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


    addpath('../Func/')
    setDir;
    load([TempDir 'TotCell.mat'], 'totCell');
    
    load([TempDir 'Paras_TWC_model.mat'], 'paras','cellToAnalysis');
    parasTWCModel     = paras;
    
    load([TempDir 'Paras_linear_model.mat'], 'paras');
    parasLinearModel  = paras;

    load([TempDir 'Paras_quadratic_model.mat'], 'paras');
    parasQuadraticModel = paras;    
    
    cellToAnalysis    = find(cellToAnalysis);
    figure;
%     
%     h1                = figure;
%     h2                = figure;

    for tCell         = 1:length(parasTWCModel)
        nCell                        = cellToAnalysis(tCell);
        para.t_frame                 = totCell(nCell).CaTime;
        para.t_ephys                 = totCell(nCell).ephysTime;
        para.fneuropil               = totCell(nCell).fNeuropil;
        para.fmean                   = totCell(nCell).fROI;
        para.filt                    = totCell(nCell).filteredEphys;
        para.peak                    = totCell(nCell).detectedSpikes;
        dff                          = get_baseline(para);
        spk                          = {para.t_ephys(para.peak)};
        para_final                   = [parasTWCModel(tCell).Fm, parasTWCModel(tCell).K, ...
                                        parasTWCModel(tCell).n, parasTWCModel(tCell).tau_d, ...
                                        parasTWCModel(tCell).tau_r];
        CaTracesTWC                  = func_getCaTraces_general_new(spk, para.t_frame, para_final);
        
        spk                          = spk{1};
        para_final                   = [parasLinearModel(tCell).a, parasLinearModel(tCell).b, ...
                                        parasLinearModel(tCell).tau_d, ...
                                        parasLinearModel(tCell).tau_r];
        CaTracesLinear               = func_getCaTraces_linear(spk, para.t_frame, para_final);
 
        para_final                   = [parasQuadraticModel(tCell).a, parasQuadraticModel(tCell).b, ...
                                        parasQuadraticModel(tCell).c, parasQuadraticModel(tCell).tau_d, ...
                                        parasQuadraticModel(tCell).tau_r];
        CaTracesQuadratic            = func_getCaTraces_quadratic(spk, para.t_frame, para_final);
        
        if strcmp(parasTWCModel(tCell).CaIndicator,'GCaMP6f')
            subplot(2,1,1)
        else
            subplot(2,1,2)
        end
        hold on;
        plot(dff, CaTracesTWC, '.k');
    end
    
    subplot(2,1,1)
    hold on;
    plot([0 14],[0 14],'--r');
    xlim([-2 14])
    ylim([0 12])
    xlabel('Baseline-corrected F')
    ylabel('F from model')
    title('GCaMP6f')
    box off
    
    subplot(2,1,2)
    hold on;
        plot([0 14],[0 14],'--r');
    xlim([-2 14])
    ylim([0 12])
    xlabel('Baseline-corrected F')
    ylabel('F from model')
    title('GCaMP6s')
    box off
    
    setPrint(8, 12, [PlotDir 'TWCModelFit'],'tif')