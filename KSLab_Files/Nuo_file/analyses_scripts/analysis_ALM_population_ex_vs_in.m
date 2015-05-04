%% supplemental figure 1

clear all
close all

% check to see if the data has been compiled
if ~exist('ALM_compiled_all_data.mat')
    
    % if not, compile all data first
    analysis_compile_all_data
    
end


load ALM_compiled_all_data

sig_selective(isnan(sig_selective))=0;


%% cell type by waveform
% waveform
spk_amp_all = min(waveform_all,[],2);
spk_width_all = [];
norm_wave_all = [];
for i_unit = 1:size(waveform_all,1)
    
    waveform_tmp = waveform_all(i_unit,:)/norm(waveform_all(i_unit,:));
    [wave_min i_peak_min] = min(waveform_tmp);
    [wave_max i_peak_max] = max(waveform_tmp(i_peak_min:end));
    spk_width_all(i_unit,1) = i_peak_max/19531*1000;
    norm_wave_all(i_unit,:) = waveform_tmp;
    
end


figure
subplot(1,2,1); hold on
plot(norm_wave_all(cellType_all(:,1)==0,:)','color',[.6 .3 0]);
plot(norm_wave_all(cellType_all(:,1)==1,:)','color','k');
plot(norm_wave_all(cellType_all(:,1)==2,:)','color',[.8 .8 .8]);


subplot(1,2,2); hold on
[y x]=hist(spk_width_all(:,1),0:.05:1);
h = bar(x,y);
set(h,'facecolor',[.6 .3 0]);
line([.45 .45],[0 max(y)*1.2],'linestyle',':','color','k');
line([.35 .35],[0 max(y)*1.2],'linestyle',':','color','k');
[y x]=hist(spk_width_all(cellType_all(:,1)==1,1),0:.05:1);
h = bar(x,y);
set(h,'facecolor','k');
[y x]=hist(spk_width_all(cellType_all(:,1)==2,1),0:.05:1);
h = bar(x,y);
set(h,'facecolor',[.8 .8 .8]);
xlim([.1 .8]);
ylim([0 1250])
xlabel('ms')
ylabel('# of cells')




%% selectivity by cell type
% pyramidal cells
cell_type_123 = zeros(size(sig_selective,1),1);
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & ~sig_selective(:,3) & cellType_all(:,1)==1)) = 1;
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==1)) = 2;
cell_type_123(find(~(sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==1)) = 3;
cell_type_tmp = [sum(cell_type_123==1) sum(cell_type_123==2) sum(cell_type_123==3)];
for i_btstrp = 1:1000
    %% sem over mice
    mice = unique(Mice_all(:,1));
    i_sample = randsample(size(mice,1),size(mice,1),'true');
    mice_iBtstrp = mice(i_sample);
    cell_type_123_iBtstrp = [];
    for i_mice = mice_iBtstrp'
        cell_type_123_iBtstrp = cat(1, cell_type_123_iBtstrp, cell_type_123(Mice_all(:,1)==i_mice,:));
    end
    cell_type_btstrp(i_btstrp,:) = [sum(cell_type_123_iBtstrp==1) sum(cell_type_123_iBtstrp==2) sum(cell_type_123_iBtstrp==3)];
    cell_type_btstrp(i_btstrp,:) = cell_type_btstrp(i_btstrp,:)/sum(cell_type_btstrp(i_btstrp,:));
end


i_sel_cell = find(sum(sig_selective,2)>0 & cellType_all(:,1)==1);
spk_count_yes_tmp = spk_count_yes_all(i_sel_cell,4);
spk_count_no_tmp = spk_count_no_all(i_sel_cell,4);
preference = [sign(spk_count_yes_tmp(:,1)-spk_count_no_tmp(:,1))];
mice_sel = Mice_all(i_sel_cell,:);
preference_btstrp = [];
for i_btstrp = 1:1000
    %% sem over mice
    mice = unique(mice_sel(:,1));
    i_sample = randsample(size(mice,1),size(mice,1),'true');
    mice_iBtstrp = mice(i_sample);
    preference_tmp = [];
    for i_mice = mice_iBtstrp'
        preference_tmp = cat(1, preference_tmp, preference(mice_sel(:,1)==i_mice,:));
    end
    preference_btstrp(i_btstrp,:) = [sum(preference_tmp>0) sum(preference_tmp<0)]/length(preference_tmp);
end
preference = [sum(preference>0) sum(preference<0)];

figure
subplot(3,2,1); hold on
h = bar(1, cell_type_tmp(1)/sum(cell_type_tmp));    set(h,'facecolor','w') 
h = bar(2, cell_type_tmp(2)/sum(cell_type_tmp));    set(h,'facecolor',[.8 .8 .8])
h = bar(3, cell_type_tmp(3)/sum(cell_type_tmp));    set(h,'facecolor','k') 
errorbar(1, cell_type_tmp(1)/sum(cell_type_tmp), std(cell_type_btstrp(:,1)),'k');
errorbar(2, cell_type_tmp(2)/sum(cell_type_tmp), std(cell_type_btstrp(:,2)),'k');
errorbar(3, cell_type_tmp(3)/sum(cell_type_tmp), std(cell_type_btstrp(:,3)),'k');
ylim([0 .6])
title(['putative pyramidal neuron n=', num2str(length(i_sel_cell))]);
legend('preparatory','prep.+ peri','peri-movement')
ylabel('Fraction of neurons')

subplot(3,2,2); hold on
bar(5, preference(1)/sum(preference),'b');
bar(6, preference(2)/sum(preference),'r');
errorbar(5, preference(1)/sum(preference), std(preference_btstrp(:,1)),'k');
errorbar(6, preference(2)/sum(preference), std(preference_btstrp(:,2)),'k');
ylim([0 .8])
legend('contra-preferring','ipsi-preferring')






% FS cells
cell_type_123 = zeros(size(sig_selective,1),1);
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & ~sig_selective(:,3) & cellType_all(:,1)==2)) = 1;
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==2)) = 2;
cell_type_123(find(~(sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==2)) = 3;
cell_type_tmp = [sum(cell_type_123==1) sum(cell_type_123==2) sum(cell_type_123==3)];
for i_btstrp = 1:1000
    %% sem over mice
    mice = unique(Mice_all(:,1));
    i_sample = randsample(size(mice,1),size(mice,1),'true');
    mice_iBtstrp = mice(i_sample);
    cell_type_123_iBtstrp = [];
    for i_mice = mice_iBtstrp'
        cell_type_123_iBtstrp = cat(1, cell_type_123_iBtstrp, cell_type_123(Mice_all(:,1)==i_mice,:));
    end
    cell_type_btstrp(i_btstrp,:) = [sum(cell_type_123_iBtstrp==1) sum(cell_type_123_iBtstrp==2) sum(cell_type_123_iBtstrp==3)];
    cell_type_btstrp(i_btstrp,:) = cell_type_btstrp(i_btstrp,:)/sum(cell_type_btstrp(i_btstrp,:));
end

i_sel_cell = find(sum(sig_selective,2)>0 & cellType_all(:,1)==2);
spk_count_yes_tmp = spk_count_yes_all(i_sel_cell,4);
spk_count_no_tmp = spk_count_no_all(i_sel_cell,4);
preference = [sign(spk_count_yes_tmp(:,1)-spk_count_no_tmp(:,1))];
mice_sel = Mice_all(i_sel_cell,:);
preference_btstrp = [];
for i_btstrp = 1:1000
    %% sem over mice
    mice = unique(mice_sel(:,1));
    i_sample = randsample(size(mice,1),size(mice,1),'true');
    mice_iBtstrp = mice(i_sample);
    preference_tmp = [];
    for i_mice = mice_iBtstrp'
        preference_tmp = cat(1, preference_tmp, preference(mice_sel(:,1)==i_mice,:));
    end
    preference_btstrp(i_btstrp,:) = [sum(preference_tmp>0) sum(preference_tmp<0)]/length(preference_tmp);
end
preference = [sum(preference>0) sum(preference<0)];

subplot(3,2,3); hold on
h = bar(1, cell_type_tmp(1)/sum(cell_type_tmp));    set(h,'facecolor','w') 
h = bar(2, cell_type_tmp(2)/sum(cell_type_tmp));    set(h,'facecolor',[.8 .8 .8])
h = bar(3, cell_type_tmp(3)/sum(cell_type_tmp));    set(h,'facecolor','k') 
errorbar(1, cell_type_tmp(1)/sum(cell_type_tmp), std(cell_type_btstrp(:,1)),'k');
errorbar(2, cell_type_tmp(2)/sum(cell_type_tmp), std(cell_type_btstrp(:,2)),'k');
errorbar(3, cell_type_tmp(3)/sum(cell_type_tmp), std(cell_type_btstrp(:,3)),'k');
ylim([0 .6])
title(['putative interneuron', num2str(length(i_sel_cell))]);

subplot(3,2,4); hold on
bar(5, preference(1)/sum(preference),'b');
bar(6, preference(2)/sum(preference),'r');
errorbar(5, preference(1)/sum(preference), std(preference_btstrp(:,1)),'k');
errorbar(6, preference(2)/sum(preference), std(preference_btstrp(:,2)),'k');
ylim([0 .8])



% unclassified cells
cell_type_123 = zeros(size(sig_selective,1),1);
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & ~sig_selective(:,3) & cellType_all(:,1)==0)) = 1;
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==0)) = 2;
cell_type_123(find(~(sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==0)) = 3;
cell_type_tmp = [sum(cell_type_123==1) sum(cell_type_123==2) sum(cell_type_123==3)];
for i_btstrp = 1:1000
    %% sem over mice
    mice = unique(Mice_all(:,1));
    i_sample = randsample(size(mice,1),size(mice,1),'true');
    mice_iBtstrp = mice(i_sample);
    cell_type_123_iBtstrp = [];
    for i_mice = mice_iBtstrp'
        cell_type_123_iBtstrp = cat(1, cell_type_123_iBtstrp, cell_type_123(Mice_all(:,1)==i_mice,:));
    end
    cell_type_btstrp(i_btstrp,:) = [sum(cell_type_123_iBtstrp==1) sum(cell_type_123_iBtstrp==2) sum(cell_type_123_iBtstrp==3)];
    cell_type_btstrp(i_btstrp,:) = cell_type_btstrp(i_btstrp,:)/sum(cell_type_btstrp(i_btstrp,:));
end

i_sel_cell = find(sum(sig_selective,2)>0 & cellType_all(:,1)==0);
spk_count_yes_tmp = spk_count_yes_all(i_sel_cell,4);
spk_count_no_tmp = spk_count_no_all(i_sel_cell,4);
preference = [sign(spk_count_yes_tmp(:,1)-spk_count_no_tmp(:,1))];
mice_sel = Mice_all(i_sel_cell,:);
preference_btstrp = [];
for i_btstrp = 1:1000
    %% sem over mice
    mice = unique(mice_sel(:,1));
    i_sample = randsample(size(mice,1),size(mice,1),'true');
    mice_iBtstrp = mice(i_sample);
    preference_tmp = [];
    for i_mice = mice_iBtstrp'
        preference_tmp = cat(1, preference_tmp, preference(mice_sel(:,1)==i_mice,:));
    end
    preference_btstrp(i_btstrp,:) = [sum(preference_tmp>0) sum(preference_tmp<0)]/length(preference_tmp);
end
preference = [sum(preference>0) sum(preference<0)];

subplot(3,2,5); hold on
h = bar(1, cell_type_tmp(1)/sum(cell_type_tmp));    set(h,'facecolor','w') 
h = bar(2, cell_type_tmp(2)/sum(cell_type_tmp));    set(h,'facecolor',[.8 .8 .8])
h = bar(3, cell_type_tmp(3)/sum(cell_type_tmp));    set(h,'facecolor','k') 
errorbar(1, cell_type_tmp(1)/sum(cell_type_tmp), std(cell_type_btstrp(:,1)),'k');
errorbar(2, cell_type_tmp(2)/sum(cell_type_tmp), std(cell_type_btstrp(:,2)),'k');
errorbar(3, cell_type_tmp(3)/sum(cell_type_tmp), std(cell_type_btstrp(:,3)),'k');
ylim([0 .6])
title(['unclassified cell', num2str(length(i_sel_cell))]);

subplot(3,2,6); hold on
bar(5, preference(1)/sum(preference),'b');
bar(6, preference(2)/sum(preference),'r');
errorbar(5, preference(1)/sum(preference), std(preference_btstrp(:,1)),'k');
errorbar(6, preference(2)/sum(preference), std(preference_btstrp(:,2)),'k');
ylim([0 .8])






%% population selectivity time course

i_sel_cell_ex = find(cellType_all(:,1)==1);
i_sel_cell_in = find(cellType_all(:,1)==2);

i_selective = sig_selective(:,1)|sig_selective(:,2)|sig_selective(:,3);
i_n_trial = N_trials_all(:,1)>15 & N_trials_all(:,2)>15;


spk_count_yes_tmp = spk_count_yes_screen(i_selective & i_n_trial,:);
spk_count_no_tmp = spk_count_no_screen(i_selective & i_n_trial,:);
preference = sign(spk_count_yes_tmp(:,1)-spk_count_no_tmp(:,1));
PSTH_prefer_tmp = PSTH_prefer_cue_aligned(i_selective & i_n_trial,:);
PSTH_nonprefer_tmp = PSTH_nonprefer_cue_aligned(i_selective & i_n_trial,:);

cellType_tmp = cellType_all(i_selective & i_n_trial,1);


figure
subplot(2,2,1); hold on
func_plot_mean_and_sem(t,PSTH_prefer_tmp(preference>0 & cellType_tmp==1,:), 'k', [.7 .7 1], 'n')
func_plot_mean_and_sem(t,PSTH_nonprefer_tmp(preference>0 & cellType_tmp==1,:), 'k', [1 .7 .7], 'n')
line([-1.3 -1.3],[0 10],'color','k')
line([-1.3 -1.3]*2,[0 10],'color','k')
line([0 0],[0 10],'color','k')
xlim([-3 2])
ylim([2 10])
xlabel('time (s)')
ylabel('selectivity (spk/s)')
title('putative pyramidal neurons')

subplot(2,2,2); hold on
func_plot_mean_and_sem(t,PSTH_nonprefer_tmp(preference<0 & cellType_tmp==1,:), 'k', [.7 .7 1], 'n')
func_plot_mean_and_sem(t,PSTH_prefer_tmp(preference<0 & cellType_tmp==1,:), 'k', [1 .7 .7], 'n')
line([-1.3 -1.3],[0 10],'color','k')
line([-1.3 -1.3]*2,[0 10],'color','k')
line([0 0],[0 10],'color','k')
xlim([-3 2])
ylim([2 10])
xlabel('time (s)')
ylabel('selectivity (spk/s)')


subplot(2,2,3); hold on
func_plot_mean_and_sem(t,PSTH_prefer_tmp(preference>0 & cellType_tmp==2,:), 'k', [.7 .7 1], 'n')
func_plot_mean_and_sem(t,PSTH_nonprefer_tmp(preference>0 & cellType_tmp==2,:), 'k', [1 .7 .7], 'n')
line([-1.3 -1.3],[0 25],'color','k')
line([-1.3 -1.3]*2,[0 25],'color','k')
line([0 0],[0 25],'color','k')
xlim([-3 2])
ylim([0 25])
xlabel('time (s)')
ylabel('selectivity (spk/s)')
title('putative FS neurons')

subplot(2,2,4); hold on
func_plot_mean_and_sem(t,PSTH_nonprefer_tmp(preference<0 & cellType_tmp==2,:), 'k', [.7 .7 1], 'n')
func_plot_mean_and_sem(t,PSTH_prefer_tmp(preference<0 & cellType_tmp==2,:), 'k', [1 .7 .7], 'n')
line([-1.3 -1.3],[0 25],'color','k')
line([-1.3 -1.3]*2,[0 25],'color','k')
line([0 0],[0 25],'color','k')
xlim([-3 2])
ylim([0 25])
xlabel('time (s)')
ylabel('selectivity (spk/s)')


