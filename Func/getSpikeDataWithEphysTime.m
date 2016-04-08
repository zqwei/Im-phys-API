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



function SpikeDataSet = getSpikeDataWithEphysTime(SpikingDataDir, SpikeFileList, minNumTrialToAnalysis, TimeToAnalysis, binSize)
    
    SpikeDataSet       = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),1000, 1);
    tot_Unit           = 0;
    
    h                  = waitbar(0,'Initializing data loads...');
    
    for nfile = 1:length(SpikeFileList)
        
        fname               = SpikeFileList(nfile).name;
        load([SpikingDataDir fname])
        
        valid_trials = (~behavior_early_report) & (behavior_report== 1) & ...
                        (task_stimulation(:,1)==0);  %#ok<NODEF>
                    
        valid_errors = (~behavior_early_report) & (behavior_report== 0) & ...
                        (task_stimulation(:,1)==0);  
                    
        numUnits     = length(neuron_single_units); %#ok<USENS>
        % numTrials    = length(behavior_report);
        task_type    = task_trial_type =='y';
        
        
        for nUnit           = 1:numUnits
            n_missing_dat   = cellfun(@isempty, neuron_single_units{nUnit});
            n_valid_yes     = valid_trials & ~n_missing_dat & task_type;
            n_valid_no      = valid_trials & ~n_missing_dat & ~task_type;
            sum_valid_yes   = sum(n_valid_yes);
            sum_valid_no    = sum(n_valid_no);

            if sum_valid_yes> minNumTrialToAnalysis && sum_valid_no> minNumTrialToAnalysis
                unit_yes_trial = nan(sum_valid_yes, length(TimeToAnalysis));
                unit_no_trial  = nan(sum_valid_no, length(TimeToAnalysis));
                unit_yes_trial_spk_time = cell(sum_valid_yes,1);
                unit_no_trial_spk_time  = cell(sum_valid_no,1);
                
                find_n_valid_yes    = find(n_valid_yes)';
                for n_trial  = 1: sum_valid_yes
                    trial_no           = find_n_valid_yes(n_trial);
                    end_delay          = task_cue_time(trial_no);
                    unit_trial_spikes  = neuron_single_units{nUnit}{trial_no};
                    unit_trial         = unit_trial_spikes - end_delay;
                    unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                         & unit_trial < max(TimeToAnalysis));
                    unit_yes_trial(n_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize;
                    unit_yes_trial_spk_time{n_trial} = unit_trial_spikes - end_delay;
                end

                find_n_valid_no    = find(n_valid_no)';

                for n_trial  = 1: sum_valid_no    
                    trial_no           = find_n_valid_no(n_trial);
                    end_delay          = task_cue_time(trial_no);
                    unit_trial_spikes  = neuron_single_units{nUnit}{trial_no};
                    unit_trial         = unit_trial_spikes - end_delay;
                    unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                         & unit_trial < max(TimeToAnalysis));
                    unit_no_trial(n_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize;  
                    unit_no_trial_spk_time{n_trial} = unit_trial_spikes - end_delay;
                end
                tot_Unit                                    = tot_Unit + 1;
                SpikeDataSet(tot_Unit).sessionIndex         = nfile;
                SpikeDataSet(tot_Unit).nUnit                = nUnit;
                SpikeDataSet(tot_Unit).unit_yes_trial       = unit_yes_trial;
                SpikeDataSet(tot_Unit).unit_yes_trial_index = find_n_valid_yes;
                SpikeDataSet(tot_Unit).unit_yes_trial_spk_time = unit_yes_trial_spk_time;
                SpikeDataSet(tot_Unit).unit_no_trial        = unit_no_trial;   
                SpikeDataSet(tot_Unit).unit_no_trial_index  = find_n_valid_no;
                SpikeDataSet(tot_Unit).unit_no_trial_spk_time  = unit_no_trial_spk_time;
                
                n_valid_yes         = valid_errors & ~n_missing_dat & task_type;
                n_valid_no          = valid_errors & ~n_missing_dat & ~task_type;
                sum_valid_yes       = sum(n_valid_yes);
                sum_valid_no        = sum(n_valid_no);
                unit_yes_error      = nan(sum_valid_yes, length(TimeToAnalysis));
                unit_no_error       = nan(sum_valid_no, length(TimeToAnalysis));
                unit_yes_error_spk_time = cell(sum_valid_yes,1);
                unit_no_error_spk_time  = cell(sum_valid_no,1);
                
                find_n_valid_yes    = find(n_valid_yes)';
                for n_trial  = 1: sum_valid_yes
                    trial_no           = find_n_valid_yes(n_trial);
                    end_delay          = task_cue_time(trial_no);
                    unit_trial_spikes  = neuron_single_units{nUnit}{trial_no};
                    unit_trial         = unit_trial_spikes - end_delay;
                    unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                         & unit_trial < max(TimeToAnalysis));
                    unit_yes_error(n_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize;
                    unit_yes_error_spk_time{n_trial} = unit_trial_spikes - end_delay;
                end

                find_n_valid_no    = find(n_valid_no)';

                for n_trial  = 1: sum_valid_no    
                    trial_no           = find_n_valid_no(n_trial);
                    end_delay          = task_cue_time(trial_no);
                    unit_trial_spikes  = neuron_single_units{nUnit}{trial_no};
                    unit_trial         = unit_trial_spikes - end_delay;
                    unit_trial         = unit_trial(unit_trial>min(TimeToAnalysis) ...
                                         & unit_trial < max(TimeToAnalysis));
                    unit_no_error(n_trial,:) = hist(unit_trial,TimeToAnalysis)/binSize; 
                    unit_no_error_spk_time{n_trial}  = unit_trial_spikes - end_delay;
                end
                SpikeDataSet(tot_Unit).unit_yes_error       = unit_yes_error;
                SpikeDataSet(tot_Unit).unit_yes_error_index = find_n_valid_yes;
                SpikeDataSet(tot_Unit).unit_yes_error_spk_time = unit_yes_error_spk_time;
                SpikeDataSet(tot_Unit).unit_no_error        = unit_no_error;   
                SpikeDataSet(tot_Unit).unit_no_error_index  = find_n_valid_no;
                SpikeDataSet(tot_Unit).unit_no_error_spk_time  = unit_no_error_spk_time;
                
                SpikeDataSet(tot_Unit).depth_in_um          = neuron_unit_info{nUnit}.depth_in_um; %#ok<USENS>
                SpikeDataSet(tot_Unit).AP_in_um             = neuron_unit_info{nUnit}.AP_ML_in_um(1);
                SpikeDataSet(tot_Unit).ML_in_um             = neuron_unit_info{nUnit}.AP_ML_in_um(2);
                SpikeDataSet(tot_Unit).cell_type            = neuron_unit_info{nUnit}.cell_type1;
                
            end
        end
        
        waitbar(nfile/length(SpikeFileList), h, sprintf('%d of %d files have been finished...',nfile, length(SpikeFileList)));
        
    end
    
    if tot_Unit < 1000
        SpikeDataSet = SpikeDataSet(1:tot_Unit);
    end
    
    close (h)