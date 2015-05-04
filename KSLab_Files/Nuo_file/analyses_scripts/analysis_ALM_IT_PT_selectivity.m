%% figure 3


clear all
close all

% check to see if the data has been compiled
if ~exist('ALM_compiled_all_data.mat')
    
    % if not, compile all data first
    analysis_compile_all_data
    
end


load ALM_compiled_all_data

sig_selective(isnan(sig_selective))=0;



%% PT neurons, plot population contra selectivity
figure; hold on
% all ALM prymidal cells
i_select_cell = (sig_selective(:,1) | sig_selective(:,2) | sig_selective(:,3)) & cellType_all(:,1)==1 & N_trials_all(:,1)>15 & N_trials_all(:,2)>15;
func_plot_mean_and_sem(t, PSTH_yes_cue_aligned(i_select_cell,:)-PSTH_no_cue_aligned(i_select_cell,:), 'k', [.8 .8 .8], 'n',3)

% PT neurons
func_plot_mean_and_sem(t, PSTH_yes_cue_aligned(cellType_all(:,2)==1,:)-PSTH_no_cue_aligned(cellType_all(:,2)==1,:), [0 .65  0], [.8 .9 .8], 'n',3)
line([0 0],[-2.5 4.5],'color','k','linestyle',':')
line([0 0]-1.3,[-2.5 4.5],'color','k','linestyle',':','color',[.6 .6 .6])
line([0 0]-2.6,[-2.5 4.5],'color','k','linestyle',':','color',[.6 .6 .6])
line([-3 2],[0 0],'color','k')
xlim([-3 2])
xlabel('Time (s)')
ylabel('Selectivity Contra-Ipsi (spk/s)')
ylim([-2.5 4.5])



%% IT neurons, plot population contra selectivity
figure; hold on
% all ALM prymidal cells
i_select_cell = (sig_selective(:,1) | sig_selective(:,2) | sig_selective(:,3)) & cellType_all(:,1)==1 & N_trials_all(:,1)>15 & N_trials_all(:,2)>15;
func_plot_mean_and_sem(t, PSTH_yes_cue_aligned(i_select_cell,:)-PSTH_no_cue_aligned(i_select_cell,:), 'k', [.8 .8 .8], 'n',3)

% IT neurons
func_plot_mean_and_sem(t, PSTH_yes_cue_aligned(cellType_all(:,2)==3,:)-PSTH_no_cue_aligned(cellType_all(:,2)==3,:), [.45 .3 .45], [.95 .85 .95], 'n',3)
line([0 0],[-2.5 4.5],'color','k','linestyle',':')
line([0 0]-1.3,[-2.5 4.5],'color','k','linestyle',':','color',[.6 .6 .6])
line([0 0]-2.6,[-2.5 4.5],'color','k','linestyle',':','color',[.6 .6 .6])
line([-3 2],[0 0],'color','k')
xlim([-3 2])
xlabel('Time (s)')
ylabel('Selectivity Contra-Ipsi (spk/s)')
ylim([-2.5 4.5])








%% fraction of neurons w/ preparatory vs. peri-movement activity
% PT cells
cell_type_123 = zeros(size(sig_selective,1),1);
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & ~sig_selective(:,3) & cellType_all(:,2)==1)) = 1;
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,2)==1)) = 2;
cell_type_123(find(~(sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,2)==1)) = 3;
cell_type_tmp = [sum(cell_type_123==1) sum(cell_type_123==2) sum(cell_type_123==3)];
cell_type_btstrp = [];
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

figure
subplot(1,2,1); hold on
h = bar(1, cell_type_tmp(1)/sum(cell_type_tmp));    set(h,'facecolor','w') 
h = bar(2, cell_type_tmp(2)/sum(cell_type_tmp));    set(h,'facecolor',[.8 .8 .8])
h = bar(3, cell_type_tmp(3)/sum(cell_type_tmp));    set(h,'facecolor','k') 
errorbar(1, cell_type_tmp(1)/sum(cell_type_tmp), std(cell_type_btstrp(:,1)),'k');
errorbar(2, cell_type_tmp(2)/sum(cell_type_tmp), std(cell_type_btstrp(:,2)),'k');
errorbar(3, cell_type_tmp(3)/sum(cell_type_tmp), std(cell_type_btstrp(:,3)),'k');
ylim([0 .8])
title('PT neurons')
legend('preparatory','prep.+ peri','peri-movement')





% IT cells
cell_type_123 = zeros(size(sig_selective,1),1);
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & ~sig_selective(:,3) & cellType_all(:,2)==3)) = 1;
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,2)==3)) = 2;
cell_type_123(find(~(sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,2)==3)) = 3;
cell_type_tmp = [sum(cell_type_123==1) sum(cell_type_123==2) sum(cell_type_123==3)];
cell_type_btstrp = [];
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

subplot(1,2,2); hold on
h = bar(1, cell_type_tmp(1)/sum(cell_type_tmp));    set(h,'facecolor','w') 
h = bar(2, cell_type_tmp(2)/sum(cell_type_tmp));    set(h,'facecolor',[.8 .8 .8])
h = bar(3, cell_type_tmp(3)/sum(cell_type_tmp));    set(h,'facecolor','k') 
errorbar(1, cell_type_tmp(1)/sum(cell_type_tmp), std(cell_type_btstrp(:,1)),'k');
errorbar(2, cell_type_tmp(2)/sum(cell_type_tmp), std(cell_type_btstrp(:,2)),'k');
errorbar(3, cell_type_tmp(3)/sum(cell_type_tmp), std(cell_type_btstrp(:,3)),'k');
ylim([0 .8])
title('IT neurons')





%% fraction of contra-preferring neurons
% all prymidal cells
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

title('All prymidal neurons')
figure; hold on
bar(5, preference(1)/sum(preference),'k');
errorbar(5, preference(1)/sum(preference), std(preference_btstrp(:,1)),'k');



% PT cells
i_sel_cell = find(sum(sig_selective,2)>0 & cellType_all(:,2)==1);
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

bar(6, preference(1)/sum(preference),'g');
errorbar(6, preference(1)/sum(preference), std(preference_btstrp(:,1)),'k');



% IT cells
i_sel_cell = find(sum(sig_selective,2)>0 & cellType_all(:,2)==3);
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

bar(7, preference(1)/sum(preference),'m');
errorbar(7, preference(1)/sum(preference), std(preference_btstrp(:,1)),'k');
ylim([0 .8])
line([4 8],[.5 .5],'linestyle',':','color','k')
ylabel('Fraction of contra-preferring neurons');
legend('all pyramidal neurons','','PT neurons','','IT neurons')



