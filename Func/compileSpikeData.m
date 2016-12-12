% 
% compile the spike raw imaging data from a list of files
% 
% version 1.0
%
% Comparison list
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function spikeDataSet    = compileSpikeData(SpikingDataDir, SpikeFileList)

    spikeDataSet         = repmat(struct('sessionName',1, 'behavior_report', 1, ...
                                'neuron_unit_info', 1, 'task_cue_time', 1),length(SpikeFileList), 1);   
    for nfile = 1:length(SpikeFileList)
        
        fname               = SpikeFileList(nfile).name;
        load([SpikingDataDir fname])
        valid_trial         = (~behavior_early_report) & ...
                              (task_stimulation(:,1)==0) & ...
                              (behavior_report>=0); %#ok<NODEF>
        spikeDataSet(nfile).sessionName      = SpikeFileList(nfile).name;
        spikeDataSet(nfile).behavior_report  = behavior_report(valid_trial);
        spikeDataSet(nfile).neuron_unit_info = neuron_unit_info;
        spikeDataSet(nfile).task_cue_time    = task_cue_time(valid_trial);
        spikeDataSet(nfile).task_pole_time   = task_pole_time(valid_trial,:);
        spikeDataSet(nfile).task_trial_type  = task_trial_type(valid_trial);
        for nNeuron         = 1:length(neuron_single_units) %#ok<USENS>
            spikeDataSet(nfile).neuron_single_units{nNeuron}  = neuron_single_units{nNeuron}(valid_trial);
        end
    end