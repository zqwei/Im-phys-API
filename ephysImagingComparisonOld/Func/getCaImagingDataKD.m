% 
% obtain the ca++ imaging data from a list of files from Kayvon's data
% 
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons) 
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function CaImagingDataSet    = getCaImagingDataKD(CaImagingDir, CaImagingFileList, minNumTrialToAnalysis, paramsROI)

    CaImagingDataSet   = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),1000, 1);
    tot_Unit           = 0;
    
    h                  = waitbar(0,'Initializing data loads...');
    
    for nfile = 1:length(CaImagingFileList)
        fname               = CaImagingFileList(nfile).name;
        load([CaImagingDir fname])
        [dFF, mean_F, correctRightTrial, correctLeftTrial, errorRightTrial, errorLeftTrial, cell_position] = ...
            extractALMImagingDataKD (dat_small, paramsROI);
        
        if sum(abs(mean_F)) == 0
            disp(['error roi file -- ' fname])
        end
        
        if length(correctRightTrial) ~= size(dFF, 3)
            disp(fname)
            continue;
        end
        
        if sum(correctRightTrial) > minNumTrialToAnalysis && sum(correctLeftTrial) > minNumTrialToAnalysis
            numUnit         = size(dFF, 1);
            dFF             = double(dFF);
            for nUnit       = 1: numUnit
                tot_Unit    = tot_Unit + 1;
                CaImagingDataSet(tot_Unit).sessionIndex         = nfile;
                CaImagingDataSet(tot_Unit).nUnit                = nUnit;
                CaImagingDataSet(tot_Unit).unit_yes_trial       = squeeze(dFF(nUnit,:,correctRightTrial))';
                CaImagingDataSet(tot_Unit).unit_yes_trial_index = find(correctRightTrial);
                CaImagingDataSet(tot_Unit).unit_no_trial        = squeeze(dFF(nUnit,:,correctLeftTrial))';
                CaImagingDataSet(tot_Unit).unit_no_trial_index  = find(correctLeftTrial);
                
                CaImagingDataSet(tot_Unit).unit_yes_error       = squeeze(dFF(nUnit,:,errorRightTrial))';
                CaImagingDataSet(tot_Unit).unit_yes_error_index = find(errorRightTrial);
                if sum(errorRightTrial) == 1
                    CaImagingDataSet(tot_Unit).unit_yes_error   = CaImagingDataSet(tot_Unit).unit_yes_error';
                end
                CaImagingDataSet(tot_Unit).unit_no_error        = squeeze(dFF(nUnit,:,errorLeftTrial))';   
                CaImagingDataSet(tot_Unit).unit_no_error_index  = find(errorLeftTrial);
                if sum(errorLeftTrial) == 1
                    CaImagingDataSet(tot_Unit).unit_no_error    = CaImagingDataSet(tot_Unit).unit_no_error';
                end
                CaImagingDataSet(tot_Unit).baseline             = mean_F(nUnit);
                
                CaImagingDataSet(tot_Unit).depth_in_um          = str2double(fname(26:28));                
  
%                 RoiCenterPos                                = ROI_list(nUnit).centerPos - (1-orgFOV)*[256 256];
%                 position_from_bregma                        = get_cell_position_from_bregma(RoiCenterPos, [ML, AP], Angle);
%                 CaImagingDataSet(tot_Unit).ML_in_um         = position_from_bregma(1);
%                 CaImagingDataSet(tot_Unit).AP_in_um         = position_from_bregma(2);
                CaImagingDataSet(tot_Unit).ML_in_um         = nan;
                CaImagingDataSet(tot_Unit).AP_in_um         = nan;
                CaImagingDataSet(tot_Unit).cell_type        = 'non_classified';
            end            
        end
        
        waitbar(nfile/length(CaImagingFileList), h, sprintf('%d of %d files have been finished...',nfile, length(CaImagingFileList)));
    end
    
    if tot_Unit < 1000
        CaImagingDataSet = CaImagingDataSet(1:tot_Unit);
    end
    
    close (h)