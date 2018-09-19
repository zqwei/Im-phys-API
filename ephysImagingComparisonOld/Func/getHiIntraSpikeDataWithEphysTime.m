% 
% obtain the spike dataset from a list of files
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



function SpikeDataSet = getHiIntraSpikeDataWithEphysTime(SpikingDataDir, SpikeFileList, TimeToAnalysis, binSize)
    
    SpikeDataSet       = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(SpikeFileList), 1);
    h                  = waitbar(0,'Initializing data loads...');
    
    
    for nfile = 1:length(SpikeFileList)        
        fname               = SpikeFileList(nfile).name;
        load([SpikingDataDir fname]);
        load([SpikingDataDir 'meta' fname]);
        
        behavior_early_report = nb.trialTypeMat(5, :) == 1;
        task_stimulation      = nb.trialTypeMat(7, :) == 1;
        valid_yes_trial       = nb.trialTypeMat(1, :) == 1 & ~behavior_early_report & ~task_stimulation;
        valid_no_trial        = nb.trialTypeMat(2, :) == 1 & ~behavior_early_report & ~task_stimulation;
        
        trialIds              = nb.trialIds;
        
        % corrects
        find_n_valid_yes      = find(valid_yes_trial)';
        
        for m_trial            = 1:sum(valid_yes_trial)
            n_trial            = find_n_valid_yes(m_trial);
            trial_no           = trialIds(n_trial);
            trial_start        = nb.trialStartTimes(n_trial);
            end_delay          = nb.trialPropertiesHash.value{3}(n_trial);
            trialTimeTag       = nb.timeSeriesArrayHash.value{1}.trial == trial_no;
            trialTimeset       = nb.timeSeriesArrayHash.value{1}.time(trialTimeTag) - trial_start;
            voltage_diff       = nb.timeSeriesArrayHash.value{1}.valueMatrix(1, trialTimeTag) ...
                                 - nb.timeSeriesArrayHash.value{1}.valueMatrix(2, trialTimeTag);
            voltage_diff(voltage_diff<20) = 0;
            [~, unit_trial_spikes] = findpeaks(voltage_diff, trialTimeset);
            
            unit_trial         = unit_trial_spikes - end_delay;
            unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                 & unit_trial < max(TimeToAnalysis));
            unit_yes_trial(m_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize;
            unit_yes_trial_spk_time{m_trial} = unit_trial_spikes - end_delay;
        end
        
        find_n_valid_no        = find(valid_no_trial)';

        for m_trial            = 1:sum(valid_no_trial)  
            n_trial            = find_n_valid_no(m_trial);
            trial_no           = trialIds(n_trial);
            trial_start        = nb.trialStartTimes(n_trial);
            end_delay          = nb.trialPropertiesHash.value{3}(n_trial);
            trialTimeTag       = nb.timeSeriesArrayHash.value{1}.trial == trial_no;
            trialTimeset       = nb.timeSeriesArrayHash.value{1}.time(trialTimeTag) - trial_start;
            voltage_diff       = nb.timeSeriesArrayHash.value{1}.valueMatrix(1, trialTimeTag) ...
                                 - nb.timeSeriesArrayHash.value{1}.valueMatrix(2, trialTimeTag);
            voltage_diff(voltage_diff<20) = 0;
            [~, unit_trial_spikes] = findpeaks(voltage_diff, trialTimeset);
            
            unit_trial         = unit_trial_spikes - end_delay;
            unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                 & unit_trial < max(TimeToAnalysis));
            unit_no_trial(m_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize;  
            unit_no_trial_spk_time{m_trial} = unit_trial_spikes - end_delay;
        end
        

        SpikeDataSet(nfile).sessionIndex         = nfile;
        SpikeDataSet(nfile).nUnit                = 1;
        SpikeDataSet(nfile).unit_yes_trial       = unit_yes_trial;
        SpikeDataSet(nfile).unit_yes_trial_index = find_n_valid_yes;
        SpikeDataSet(nfile).unit_yes_trial_spk_time = unit_yes_trial_spk_time;
        SpikeDataSet(nfile).unit_no_trial        = unit_no_trial;   
        SpikeDataSet(nfile).unit_no_trial_index  = find_n_valid_no;
        SpikeDataSet(nfile).unit_no_trial_spk_time  = unit_no_trial_spk_time;

        % errors
        valid_yes_error       = nb.trialTypeMat(3, :) == 1 & ~behavior_early_report & ~task_stimulation;
        valid_no_error        = nb.trialTypeMat(4, :) == 1 & ~behavior_early_report & ~task_stimulation;
        
        find_n_valid_yes      = find(valid_yes_error)';
        
        for m_trial            = 1:sum(valid_yes_error)
            n_trial            = find_n_valid_yes(m_trial);
            trial_no           = trialIds(n_trial);
            trial_start        = nb.trialStartTimes(n_trial);
            end_delay          = nb.trialPropertiesHash.value{3}(n_trial);
            trialTimeTag       = nb.timeSeriesArrayHash.value{1}.trial == trial_no;
            trialTimeset       = nb.timeSeriesArrayHash.value{1}.time(trialTimeTag) - trial_start;
            voltage_diff       = nb.timeSeriesArrayHash.value{1}.valueMatrix(1, trialTimeTag) ...
                                 - nb.timeSeriesArrayHash.value{1}.valueMatrix(2, trialTimeTag);
            voltage_diff(voltage_diff<20) = 0;
            [~, unit_trial_spikes] = findpeaks(voltage_diff, trialTimeset);
            
            unit_trial         = unit_trial_spikes - end_delay;
            unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                 & unit_trial < max(TimeToAnalysis));
            unit_yes_error(m_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize;
            unit_yes_error_spk_time{m_trial} = unit_trial_spikes - end_delay;
        end
        
        find_n_valid_no        = find(valid_no_error)';

        for m_trial            = 1:sum(valid_no_error)  
            n_trial            = find_n_valid_no(m_trial);
            trial_no           = trialIds(n_trial);
            trial_start        = nb.trialStartTimes(n_trial);
            end_delay          = nb.trialPropertiesHash.value{3}(n_trial);
            trialTimeTag       = nb.timeSeriesArrayHash.value{1}.trial == trial_no;
            trialTimeset       = nb.timeSeriesArrayHash.value{1}.time(trialTimeTag) - trial_start;
            voltage_diff       = nb.timeSeriesArrayHash.value{1}.valueMatrix(1, trialTimeTag) ...
                                 - nb.timeSeriesArrayHash.value{1}.valueMatrix(2, trialTimeTag);
            voltage_diff(voltage_diff<20) = 0;
            [~, unit_trial_spikes] = findpeaks(voltage_diff, trialTimeset);
            
            unit_trial         = unit_trial_spikes - end_delay;
            unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                 & unit_trial < max(TimeToAnalysis));
            unit_no_error(m_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize;  
            unit_no_error_spk_time{m_trial} = unit_trial_spikes - end_delay;
        end
               
        
        SpikeDataSet(nfile).unit_yes_error       = unit_yes_error;
        SpikeDataSet(nfile).unit_yes_error_index = find_n_valid_yes;
        SpikeDataSet(nfile).unit_yes_error_spk_time = unit_yes_error_spk_time;
        SpikeDataSet(nfile).unit_no_error        = unit_no_error;   
        SpikeDataSet(nfile).unit_no_error_index  = find_n_valid_no;
        SpikeDataSet(nfile).unit_no_error_spk_time  = unit_no_error_spk_time;
        
        location_string                          = meta_data.intracellular.recordingCoordinates{1};
        locations                                = sscanf(location_string,'%fmm anterior, %fmm lateral, %fmm deep, ');
        SpikeDataSet(nfile).depth_in_um          = locations(3)*1000;
        SpikeDataSet(nfile).AP_in_um             = locations(1)*1000;
        SpikeDataSet(nfile).ML_in_um             = locations(2)*1000;
        SpikeDataSet(nfile).cell_type            = 2 - strcmp(meta_data.intracellular.cellType{1}, 'Excitatory');
        SpikeDataSet(nfile).animalNameIndex      = meta_data.animalID;
        
        waitbar(nfile/length(SpikeFileList), h, sprintf('%d of %d files have been finished...',nfile, length(SpikeFileList)));
        
    end
    
    
    close (h)