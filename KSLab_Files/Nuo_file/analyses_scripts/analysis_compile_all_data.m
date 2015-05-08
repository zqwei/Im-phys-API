clear all
close all

file_list = dir([file_dir,'data_structure_*.mat']);


%% establish variables to compute and save
N_trials_all = [];

spike_times_all = {};       %[yes no yes_error no_error]

PSTH_yes_cue_aligned  = [];
PSTH_no_cue_aligned  = [];
PSTH_yes_error_cue_aligned  = [];
PSTH_no_error_cue_aligned  = [];

spk_count_yes_all  = [];    %[sample delay resposne whole_trial baseline]
spk_count_no_all  = [];
spk_count_yes_error_all  = [];
spk_count_no_error_all  = [];

sig_selective = [];         %[sample delay resposne]


spk_count_yes_screen = [];
spk_count_no_screen = [];
PSTH_prefer_cue_aligned  = [];
PSTH_nonprefer_cue_aligned  = [];
Dprime_cue_aligned  = [];


Probe_type_all = {};
Channel_all = [];
cellType_all = [];
waveform_all = [];
Mice_all = [];




%% go through the files and compute the variables
n_file = length(file_list);
animal_ID = {};
n_mice = 0;
n_unit = 0;
for i_file = 1:n_file
        
    % --------------- load 'obj' and 'meta_data' -------------------
    clear obj meta_data
    dir_tmp = '..\datafiles\';
    filename_tmp = file_list(i_file).name;
    i_str = findstr(filename_tmp,'ANM');
    data_struct_filename_tmp = ['data_structure_',filename_tmp(i_str(1):end)];
    meta_data_filename_tmp = ['meta_data_',filename_tmp(i_str(1):end)];
    disp([filename_tmp(i_str(1):end),'  File #',num2str(i_file),'/',num2str(n_file)]);
    load([dir_tmp,data_struct_filename_tmp]);
    load([dir_tmp,meta_data_filename_tmp]);
            

    % get mouse ID
    animal_ID_iFile = filename_tmp(i_str(1):i_str(1)+8);
    if sum(strcmp(animal_ID,animal_ID_iFile))==0
        n_mice = n_mice+1;
        animal_ID{n_mice,1} = animal_ID_iFile;
    end
    
    
    
    % --------------- get behavioral data -----------------
    clear trialID i_Hit* i_Err* i_NoLick* i_LickEarly* i_StimTrials* i_GoodTrials
    clear Hit* Err* NoLick* LickEarly* StimTrials* GoodTrials
    
    % get the raw data
    trialID = obj.trialIds;
    
    i_HitR = find(strcmp(obj.trialTypeStr,'HitR'));
    i_HitL = find(strcmp(obj.trialTypeStr,'HitL'));
    i_ErrR = find(strcmp(obj.trialTypeStr,'ErrR'));
    i_ErrL = find(strcmp(obj.trialTypeStr,'ErrL'));
    i_NoLickR = find(strcmp(obj.trialTypeStr,'NoLickR'));
    i_NoLickL = find(strcmp(obj.trialTypeStr,'NoLickL'));
    i_LickEarly = find(strcmp(obj.trialTypeStr,'LickEarly'));
    i_StimTrials = find(strcmp(obj.trialTypeStr,'StimTrials'));
    
    HitR = obj.trialTypeMat(i_HitR,:)';
    HitL = obj.trialTypeMat(i_HitL,:)';
    ErrR = obj.trialTypeMat(i_ErrR,:)';
    ErrL = obj.trialTypeMat(i_ErrL,:)';
    NoLickR = obj.trialTypeMat(i_NoLickR,:)';
    NoLickL = obj.trialTypeMat(i_NoLickL,:)';
    LickEarly = obj.trialTypeMat(i_LickEarly,:)';
    StimTrials = obj.trialTypeMat(i_StimTrials,:)';
    
    
    % get the good trials (periods when mice is performing)
    i_GoodTrials = find(strcmp(obj.trialPropertiesHash.keyNames,'GoodTrials'));
    GoodTrials = obj.trialPropertiesHash.value{i_GoodTrials};
    
    
    % get trial temporal structure
    t_TrialStart = obj.trialStartTimes;
    
    Pole_in_Time = obj.trialPropertiesHash.value{1};
    Pole_out_Time = obj.trialPropertiesHash.value{2};
    Cue_Time = obj.trialPropertiesHash.value{3};
    
    

    % get photostimulation information 
    Photostim_type = obj.trialPropertiesHash.value{5};      % 0--no stim; 1--PT axonal stim; 2--IT axonal stim; Nan--DON'T ANALYZE (other test trials)
            
    
    % -------------- process single unit data --------------
    for i_unit = 1:length(obj.eventSeriesHash.value)
        
        unit = obj.eventSeriesHash.value{i_unit};
                
        n_unit = n_unit+1;
        
        
        % ---------- trial selection -----------------
        % index stable trials
        unit_stable_trials_tmp = zeros(trialID(end),1);
        i_stable_trials = min(unit.eventTrials):min([max(unit.eventTrials) trialID(end)]);
        unit_stable_trials_tmp(i_stable_trials) = 1;
        
        
                    
        % index the behavioral trials
        i_yes_trial_correct = find(HitR==1 & Photostim_type==0 & unit_stable_trials_tmp==1 & LickEarly==0 & GoodTrials==1);
        i_no_trial_correct = find(HitL==1 & Photostim_type==0 & unit_stable_trials_tmp==1 & LickEarly==0 & GoodTrials==1);
        i_yes_trial_error = find((ErrR==1) & Photostim_type==0 & unit_stable_trials_tmp==1 & LickEarly==0 & GoodTrials==1);
        i_no_trial_error = find((ErrL==1) & Photostim_type==0 & unit_stable_trials_tmp==1 & LickEarly==0 & GoodTrials==1);
        
        i_stim_tagging = find((Photostim_type==1 | Photostim_type==2 | Photostim_type==3 | Photostim_type==4 | Photostim_type==5) & unit_stable_trials_tmp==1 & GoodTrials==1);
        
        
        N_trials_all = [N_trials_all; length(i_yes_trial_correct) length(i_no_trial_correct) length(i_yes_trial_error) length(i_no_trial_error)  length(i_stim_tagging)];
        
        
        % ============  spike times & PSTH  ============
        % yes trials
        spike_times_psth = {};
        n_trial = 0;
        for i_trial = i_yes_trial_correct'
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);            
            spike_times_psth{n_trial,1} = spike_times_iTrial - Cue_Time(i_trial,1);
        end
        [psth1 t] = func_getPSTH(spike_times_psth,-3.5,2);
        spike_times_yes = spike_times_psth;
        
                    
        % no trials
        spike_times_psth = {};
        n_trial = 0;
        for i_trial = i_no_trial_correct'
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spike_times_psth{n_trial,1} = spike_times_iTrial - Cue_Time(i_trial,1);
        end
        [psth2 t] = func_getPSTH(spike_times_psth,-3.5,2);
        spike_times_no = spike_times_psth;
        
                    
        % yes error trials
        spike_times_psth = {};
        n_trial = 0;
        for i_trial = i_yes_trial_error'
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spike_times_psth{n_trial,1} = spike_times_iTrial - Cue_Time(i_trial,1);
        end
        [psth3 t] = func_getPSTH(spike_times_psth,-3.5,2);
        spike_times_yes_error = spike_times_psth;
                    
                    
        % no error trials
        spike_times_psth = {};
        n_trial = 0;
        for i_trial = i_no_trial_error'
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spike_times_psth{n_trial,1} = spike_times_iTrial - Cue_Time(i_trial,1);
        end
        [psth4 t] = func_getPSTH(spike_times_psth,-3.5,2);
        spike_times_no_error = spike_times_psth;
        
                    
        % save the data
        spike_times_all{n_unit,1} = spike_times_yes;
        spike_times_all{n_unit,2} = spike_times_no;
        spike_times_all{n_unit,3} = spike_times_yes_error;
        spike_times_all{n_unit,4} = spike_times_no_error;
        
                    
        PSTH_yes_cue_aligned(n_unit,:) = psth1;
        PSTH_no_cue_aligned(n_unit,:) = psth2;
        PSTH_yes_error_cue_aligned(n_unit,:) = psth3;
        PSTH_no_error_cue_aligned(n_unit,:) = psth4;
        
        
                    
                    
                    
        
        % ============ spike counts ========
        % yes trials correct
        n_trial = 0;
        spk_count_yes = [];
        for i_trial = i_yes_trial_correct'
            
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_yes(n_trial,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Pole_out_Time(i_trial,1))/(Pole_out_Time(i_trial,1)-Pole_in_Time(i_trial,1));
            spk_count_yes(n_trial,2) = sum(spike_times_iTrial>Pole_out_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1))/(Cue_Time(i_trial,1)-Pole_out_Time(i_trial,1));
            spk_count_yes(n_trial,3) = sum(spike_times_iTrial>Cue_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(1.3);
            spk_count_yes(n_trial,4) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(Cue_Time(i_trial,1)+1.3-Pole_in_Time(i_trial,1));
            spk_count_yes(n_trial,5) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1)-.5 & spike_times_iTrial<Pole_in_Time(i_trial,1))/(.5);
            
        end
                    
                    
        % no trials correct
        n_trial = 0;
        spk_count_no = [];
        for i_trial = i_no_trial_correct'
            
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_no(n_trial,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Pole_out_Time(i_trial,1))/(Pole_out_Time(i_trial,1)-Pole_in_Time(i_trial,1));
            spk_count_no(n_trial,2) = sum(spike_times_iTrial>Pole_out_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1))/(Cue_Time(i_trial,1)-Pole_out_Time(i_trial,1));
            spk_count_no(n_trial,3) = sum(spike_times_iTrial>Cue_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(1.3);
            spk_count_no(n_trial,4) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(Cue_Time(i_trial,1)+1.3-Pole_in_Time(i_trial,1));
            spk_count_no(n_trial,5) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1)-.5 & spike_times_iTrial<Pole_in_Time(i_trial,1))/(.5);
            
        end
        
                    
                    
        % yes trials error
        n_trial = 0;
        spk_count_yes_error = [];
        for i_trial = i_yes_trial_error'
            
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_yes_error(n_trial,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Pole_out_Time(i_trial,1))/(Pole_out_Time(i_trial,1)-Pole_in_Time(i_trial,1));
            spk_count_yes_error(n_trial,2) = sum(spike_times_iTrial>Pole_out_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1))/(Cue_Time(i_trial,1)-Pole_out_Time(i_trial,1));
            spk_count_yes_error(n_trial,3) = sum(spike_times_iTrial>Cue_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(1.3);
            spk_count_yes_error(n_trial,4) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(Cue_Time(i_trial,1)+1.3-Pole_in_Time(i_trial,1));
            spk_count_yes_error(n_trial,5) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1)-.5 & spike_times_iTrial<Pole_in_Time(i_trial,1))/(.5);
            
        end
                    
        
        % no trials error
        n_trial = 0;
        spk_count_no_error = [];
        for i_trial = i_no_trial_error'
            
            n_trial = n_trial+1;
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_no_error(n_trial,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Pole_out_Time(i_trial,1))/(Pole_out_Time(i_trial,1)-Pole_in_Time(i_trial,1));
            spk_count_no_error(n_trial,2) = sum(spike_times_iTrial>Pole_out_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1))/(Cue_Time(i_trial,1)-Pole_out_Time(i_trial,1));
            spk_count_no_error(n_trial,3) = sum(spike_times_iTrial>Cue_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(1.3);
            spk_count_no_error(n_trial,4) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1)+1.3)/(Cue_Time(i_trial,1)+1.3-Pole_in_Time(i_trial,1));
            spk_count_no_error(n_trial,5) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1)-.5 & spike_times_iTrial<Pole_in_Time(i_trial,1))/(.5);
            
        end
        
        
        % save the data
        spk_count_yes_all(n_unit,1:5) = mean(spk_count_yes);
        spk_count_no_all(n_unit,1:5) = mean(spk_count_no);
        spk_count_yes_error_all(n_unit,1:5) = mean(spk_count_yes_error);
        spk_count_no_error_all(n_unit,1:5) = mean(spk_count_no_error);
        
                    
                    
                    
                    
                    
                    
        % ============ significant selectivity ========
        % sample period, correct trials
        spk_count_yes = [];
        for i_trial = i_yes_trial_correct'
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_yes(end+1,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Pole_out_Time(i_trial,1))/(Pole_out_Time(i_trial,1)-Pole_in_Time(i_trial,1));
        end
        spk_count_no = [];
        for i_trial = i_no_trial_correct'
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_no(end+1,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<Pole_out_Time(i_trial,1))/(Pole_out_Time(i_trial,1)-Pole_in_Time(i_trial,1));
        end
        if ~isempty(spk_count_yes) & ~isempty(spk_count_no)
            sample_selective = ttest2(spk_count_yes,spk_count_no);
        else
            sample_selective = nan;
        end
                    
                    
        % delay period, correct trials
        spk_count_yes = [];
        for i_trial = i_yes_trial_correct'
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_yes(end+1,1) = sum(spike_times_iTrial>Pole_out_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1))/(Cue_Time(i_trial,1)-Pole_out_Time(i_trial,1));
        end
        spk_count_no = [];
        for i_trial = i_no_trial_correct'
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_no(end+1,1) = sum(spike_times_iTrial>Pole_out_Time(i_trial,1) & spike_times_iTrial<Cue_Time(i_trial,1))/(Cue_Time(i_trial,1)-Pole_out_Time(i_trial,1));
        end
        if ~isempty(spk_count_yes) & ~isempty(spk_count_no)
            delay_selective = ttest2(spk_count_yes,spk_count_no);
        else
            delay_selective = nan;
        end
                    
                    
        % response period, correct trials
        spk_count_yes = [];
        for i_trial = i_yes_trial_correct'
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_yes(end+1,1) = sum(spike_times_iTrial>Cue_Time(i_trial,1) & spike_times_iTrial<(Cue_Time(i_trial,1)+1.3))/1.3;
        end
        spk_count_no = [];
        for i_trial = i_no_trial_correct'
            spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
            spk_count_no(end+1,1) = sum(spike_times_iTrial>Cue_Time(i_trial,1) & spike_times_iTrial<(Cue_Time(i_trial,1)+1.3))/1.3;
        end
        if ~isempty(spk_count_yes) & ~isempty(spk_count_no)
            response_selective = ttest2(spk_count_yes,spk_count_no);
        else
            response_selective = nan;
        end
        
        
        
        % save the data
        sig_selective(n_unit,1:3) = [sample_selective delay_selective response_selective];
        
        
        
        
        
        
        % ============  PSTH preferred & nonpreferred  ============
        if length(i_yes_trial_correct)>15 & length(i_no_trial_correct)>15
            
            % select screen and test trials
            n_trials_tmp = length(i_yes_trial_correct);
            i_yes_trial_correct = i_yes_trial_correct(randsample(n_trials_tmp,n_trials_tmp));
            i_yes_trial_screen = i_yes_trial_correct(1:10);
            i_yes_trial_test = i_yes_trial_correct(11:end);
            
            n_trials_tmp = length(i_no_trial_correct);
            i_no_trial_correct = i_no_trial_correct(randsample(n_trials_tmp,n_trials_tmp));
            i_no_trial_screen = i_no_trial_correct(1:10);
            i_no_trial_test = i_no_trial_correct(11:end);
            
            
            % yes trials  screen
            spk_count_yes = [];
            n_trial = 0;
            for i_trial = i_yes_trial_screen'
                spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
                spk_count_yes(end+1,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<(Cue_Time(i_trial,1)+1.3));
            end
            
            
            % no trials  screen
            spk_count_no = [];
            n_trial = 0;
            for i_trial = i_no_trial_screen'
                spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
                spk_count_no(end+1,1) = sum(spike_times_iTrial>Pole_in_Time(i_trial,1) & spike_times_iTrial<(Cue_Time(i_trial,1)+1.3));
            end
            
                        
                        
                        
            % compute PSTH with remaining data
            % yes trials test
            spike_times_psth = {};
            n_trial = 0;
            for i_trial = i_yes_trial_test'
                n_trial = n_trial+1;
                spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
                spike_times_psth{n_trial,1} = spike_times_iTrial - Cue_Time(i_trial,1);
            end
            [psth1 t] = func_getPSTH(spike_times_psth,-3.5,2);
            
            
            
            % no trials test
            spike_times_psth = {};
            n_trial = 0;
            for i_trial = i_no_trial_test'
                n_trial = n_trial+1;
                spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
                spike_times_psth{n_trial,1} = spike_times_iTrial - Cue_Time(i_trial,1);
            end
            [psth2 t] = func_getPSTH(spike_times_psth,-3.5,2);
            
                        
                        
            % selectivity
            Dprime_iCell = [];
            Dprime_t = [];
            for i_time = (Cue_Time(i_trial,1)-3.5):.1:(Cue_Time(i_trial,1)+2)
                t_win_on = i_time;
                t_win_off = i_time+.1;
                spk_count_yes_iTime = [];
                for i_trial = i_yes_trial_test'
                    spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
                    spk_count_yes_iTime(end+1,1) = sum(spike_times_iTrial>t_win_on & spike_times_iTrial<t_win_off);
                end
                spk_count_no_iTime = [];
                for i_trial = i_no_trial_test'
                    spike_times_iTrial = unit.eventTimes(unit.eventTrials==i_trial) - t_TrialStart(i_trial);
                    spk_count_no_iTime(end+1,1) = sum(spike_times_iTrial>t_win_on & spike_times_iTrial<t_win_off);
                end
                Dprime_iCell(end+1,1) = (mean(spk_count_yes_iTime)-mean(spk_count_no_iTime))/sqrt(var(spk_count_yes_iTime)/2+var(spk_count_no_iTime)/2);
                Dprime_t(end+1,1) = i_time-Cue_Time(i_trial,1);
            end
            Dprime_iCell(isnan(Dprime_iCell)) = 0;
                        
                        
        else
            
            % if unit has <15 trials
            spk_count_yes = nan;
            spk_count_no = nan;
            Dprime_iCell = nan(1,56);
            psth1 = nan(1,5101);
            psth2 = nan(1,5101);
            
        end
                    
        
        
        % save data
        spk_count_yes_screen(n_unit,1) = mean(spk_count_yes);
        spk_count_no_screen(n_unit,1) = mean(spk_count_no);
        Dprime_cue_aligned(n_unit,:) = Dprime_iCell;
        if mean(spk_count_yes)>mean(spk_count_no)
            PSTH_prefer_cue_aligned(n_unit,:) = psth1;
            PSTH_nonprefer_cue_aligned(n_unit,:) = psth2;
        else
            PSTH_prefer_cue_aligned(n_unit,:) = psth2;
            PSTH_nonprefer_cue_aligned(n_unit,:) = psth1;
        end
        
                        
        
                    
        % ============ other  information ==============
        Probe_type_all{n_unit,1} = meta_data.extracellular.probeType;
        Channel_all(n_unit,1) = median(unit.channel);
        
        if sum(strcmp(unit.cellType,'pyramidal'))>0 & sum(strcmp(unit.cellType,'FS'))==0
            cellType_all(n_unit,1) = 1;
        elseif sum(strcmp(unit.cellType,'pyramidal'))==0 & sum(strcmp(unit.cellType,'FS'))>0
            cellType_all(n_unit,1) = 2;
        elseif sum(strcmp(unit.cellType,'pyramidal'))==0 & sum(strcmp(unit.cellType,'FS'))==0
            cellType_all(n_unit,1) = 0;
        else
            error('cell type could not be processed')
        end
        
        if sum(strcmp(unit.cellType,'PT'))>0 & sum(strcmp(unit.cellType,'IT'))==0
            cellType_all(n_unit,2) = 1;
        elseif sum(strcmp(unit.cellType,'PT'))==0 & sum(strcmp(unit.cellType,'IT'))>0
            cellType_all(n_unit,2) = 3;
        elseif sum(strcmp(unit.cellType,'PT'))==0 & sum(strcmp(unit.cellType,'IT'))==0
            cellType_all(n_unit,2) = 0;
        else
            error('cell type could not be processed')
        end
        
        waveform_all(n_unit,:) = mean(unit.waveforms);
        
        Mice_all(n_unit,:) = [n_mice    meta_data.extracellular.penetrationN    i_unit];
        
        
    end
end
clear obj


save ALM_compiled_all_data













