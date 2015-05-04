%% figure 2,  supplemental figure 1, supplemental figure 2


clear all
close all

% check to see if the data has been compiled
if ~exist('ALM_compiled_all_data.mat')
    
    % if not, compile all data first
    analysis_compile_all_data
    
end


load ALM_compiled_all_data

sig_selective(isnan(sig_selective))=0;



%% example neurons, PSTH and raster
example_cells = {
    'ANM218457', 4, 9;...
    'ANM219248', 2, 11;...
    'ANM219037', 2, 5;...
    'ANM219031', 2, 17;...
    'ANM210862', 2, 21;...
    'ANM210863', 2, 10;...
    'ANM210863', 1, 9;...
    'ANM219033', 5, 1;...
    'ANM219036', 5, 13;...
    };

figure
n_plot = 0;
for i_ex_cell = 1:size(example_cells,1)
    i_mice = find(strcmp(animal_ID,example_cells{i_ex_cell,1}));
    i_session = example_cells{i_ex_cell,2};
    i_unit = example_cells{i_ex_cell,3};
    n_plot = n_plot+1;
    
    i_cell = find(Mice_all(:,1)==i_mice & Mice_all(:,2)==i_session & Mice_all(:,3)==i_unit);
    spike_times1 = spike_times_all{i_cell,1};
    spike_times2 = spike_times_all{i_cell,2};
    
    psth1 = PSTH_yes_cue_aligned(i_cell,:);
    psth2 = PSTH_no_cue_aligned(i_cell,:);
    
    subplot(3,3,n_plot); hold on
    plot(t,psth2,'r');
    plot(t,psth1,'b');
    line([0 0],[0 2.2]*max([psth1 psth2]),'color','k')
    line([-1.3 -1.3],[0 2.2]*max([psth1 psth2]),'color','k')
    line([-2.6 -2.6],[0 2.2]*max([psth1 psth2]),'color','k')
    xlim([-3 2])
    
    y_max = max([psth1 psth2])*1.2;
    y_scale = max([psth1 psth2]);
    n_trials = size(spike_times1,1)+size(spike_times2,1);
    
    for i=1:length(spike_times2)
        if ~isempty(spike_times2{i})
            line([spike_times2{i} spike_times2{i}]', [y_max+y_scale/n_trials*(i-1) y_max+y_scale/n_trials*i],'color','r')            
        end
    end
    y_max = y_max+y_scale/n_trials*i;
    for i=1:length(spike_times1)
        if ~isempty(spike_times1{i})
            line([spike_times1{i} spike_times1{i}]', [y_max+y_scale/n_trials*(i-1) y_max+y_scale/n_trials*i],'color','b')
        end
    end
    
    ylim([0 2.3]*max([psth1 psth2]));
    
    if n_plot==1
        xlabel('Time (s)')
        ylabel('Spk/s')
        legend('lick left','lick right')
    end
    
end






%% example neurons, PSTH (supplemental figure 2)

example_cells = {
    'ANM214430', 2, 4;...
    'ANM218457', 1, 1;...
    'ANM210863', 1, 7;...
    'ANM219033', 4, 5;...
    'ANM210862', 2, 16;...
    'ANM210861', 2, 5;...
    'ANM210863', 1, 10;...
    'ANM218457', 3, 3;...
    'ANM219031', 5, 3;...
    'ANM210863', 2, 10;...
    'ANM219033', 3, 5;...
    'ANM218457', 4, 12;...
    'ANM219247', 1, 3;...
    'ANM219248', 2, 11;...
    'ANM210862', 2, 4;...
    'ANM219248', 1, 11;...
    'ANM219031', 5, 4;...
    'ANM219031', 2, 14;...
    'ANM210863', 3, 33;...    
    'ANM219030', 5, 8;...
    };

figure
n_plot = 0;
for i_ex_cell = 1:size(example_cells,1)
    i_mice = find(strcmp(animal_ID,example_cells{i_ex_cell,1}));
    i_session = example_cells{i_ex_cell,2};
    i_unit = example_cells{i_ex_cell,3};
    n_plot = n_plot+1;
    
    i_cell = find(Mice_all(:,1)==i_mice & Mice_all(:,2)==i_session & Mice_all(:,3)==i_unit);
    spike_times1 = spike_times_all{i_cell,1};
    spike_times2 = spike_times_all{i_cell,2};
    
    psth1 = PSTH_yes_cue_aligned(i_cell,:);
    psth2 = PSTH_no_cue_aligned(i_cell,:);
    
    subplot(4,5,n_plot); hold on
    plot(t,psth2,'r');
    plot(t,psth1,'b');
    line([0 0],[0 2.2]*max([psth1 psth2]),'color','k')
    line([-1.3 -1.3],[0 2.2]*max([psth1 psth2]),'color','k')
    line([-2.6 -2.6],[0 2.2]*max([psth1 psth2]),'color','k')
    xlim([-3 2])
    
    ylim([0 2.3]*max([psth1 psth2]));
    
    if n_plot==1
        xlabel('Time (s)')
        ylabel('Spk/s')
        legend('lick left','lick right')
    end
    
end




%% poulation selectivity time course

% sorted by contra and ipsi, prymidal cells only
FR_diff_all = PSTH_yes_cue_aligned-PSTH_no_cue_aligned;
FR_diff_all = -FR_diff_all;  % flip color so blue means contra

% preparatory, prep.+peri, peri-movement
i_type_I = (sig_selective(:,1)|sig_selective(:,2)) & ~sig_selective(:,3) & cellType_all(:,1)==1;
i_type_II = (sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==1;
i_type_III = ~(sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==1;


% sort the neurons by their preparatory activity
preference = [];
a=[];
cell_type123 = [];

% preparatory only
FR_diff_all_tmp = FR_diff_all(i_type_I,:);
spk_count_yes_tmp = spk_count_yes_all(i_type_I,4);
spk_count_no_tmp = spk_count_no_all(i_type_I,4);
for i_trial = 1:size(FR_diff_all_tmp,1)
    FR_diff_all_tmp(i_trial,:) = FR_diff_all_tmp(i_trial,:);
    FR_diff_all_tmp(i_trial,:)  = FR_diff_all_tmp(i_trial,:) / max(abs(FR_diff_all_tmp(i_trial,:)));
end
preference = [preference; sign(spk_count_yes_tmp(:,1)-spk_count_no_tmp(:,1))];
a = [a; FR_diff_all_tmp(:,:)];
cell_type123 = [cell_type123; ones(sum(i_type_I),1)]


% preparatory + perimovement
FR_diff_all_tmp = FR_diff_all(i_type_II,:);
spk_count_yes_tmp = spk_count_yes_all(i_type_II,4);
spk_count_no_tmp = spk_count_no_all(i_type_II,4);
for i_trial = 1:size(FR_diff_all_tmp,1)
    FR_diff_all_tmp(i_trial,:) = FR_diff_all_tmp(i_trial,:);
    FR_diff_all_tmp(i_trial,:)  = FR_diff_all_tmp(i_trial,:) / max(abs(FR_diff_all_tmp(i_trial,:)));
end
preference = [preference; sign(spk_count_yes_tmp(:,1)-spk_count_no_tmp(:,1))];
a = [a; FR_diff_all_tmp(:,:)];
cell_type123 = [cell_type123; ones(sum(i_type_II),1)*2]


% peri-movement
FR_diff_all_tmp = FR_diff_all(i_type_III,:);
spk_count_yes_tmp = spk_count_yes_all(i_type_III,4);
spk_count_no_tmp = spk_count_no_all(i_type_III,4);
for i_trial = 1:size(FR_diff_all_tmp,1)
    FR_diff_all_tmp(i_trial,:) = FR_diff_all_tmp(i_trial,:);
    FR_diff_all_tmp(i_trial,:)  = FR_diff_all_tmp(i_trial,:) / max(abs(FR_diff_all_tmp(i_trial,:)));
end
FR_diff_all_tmp(isnan(FR_diff_all_tmp(:,1)),:)=[];
preference = [preference; sign(spk_count_yes_tmp(:,1)-spk_count_no_tmp(:,1))];
a = [a; FR_diff_all_tmp(:,:)];
cell_type123 = [cell_type123; ones(sum(i_type_III),1)*3]



% sort the neurons by contra-preferring vs. ipsi-preferring
i_contra = find(preference>0);
i_ipsi = find(preference<0);

figure
subplot(2,10,[1:9]);
y=linspace(1,1,size(i_contra,2));
x=t;
imagesc(x,y,a(i_contra,:))
line([0 0],[1 size(a,2)],'color','k')
line([-1.3 -1.3],[1 size(a,2)],'color','k')
line([-1.3 -1.3]*2,[1 size(a,2)],'color','k')
xlabel('time (s)')
xlim([-3 1.8]) 

subplot(2,10,10);
imagesc(cell_type123(i_contra));



subplot(2,10,[11:19]);
y=linspace(1,1,size(i_ipsi,2));
x=t;
imagesc(x,y,a(i_ipsi,:))
line([0 0],[1 size(a,2)],'color','k')
line([-1.3 -1.3],[1 size(a,2)],'color','k')
line([-1.3 -1.3]*2,[1 size(a,2)],'color','k')
xlabel('time (s)')
ylabel('cell #')
xlim([-3 1.8]) 

subplot(2,10,20);
imagesc(cell_type123(i_ipsi));





%% population average selectivity time course 
i_sel = (cellType_all(:,1)==1 & sum(sig_selective,2)>0 & N_trials_all(:,1)>15 & N_trials_all(:,2)>15 & (spk_count_yes_all(:,4)>.1 & spk_count_no_all(:,4)>.1));
FR_diff_all = PSTH_prefer_cue_aligned-PSTH_nonprefer_cue_aligned;
FR_diff_all = FR_diff_all(i_sel,:);


figure
subplot(2,1,1);
func_plot_mean_and_sem(t, FR_diff_all, 'k', [.8 .8 .8], 'n')
y_max = max([mean(FR_diff_all)]);
y_min = min([mean(FR_diff_all)]);
line([0 0],[-1 4],'color','k','linestyle',':')
line([0 0]-1.3,[-1 4],'color','k','linestyle',':','color',[.6 .6 .6])
line([0 0]-2.6,[-1 4],'color','k','linestyle',':','color',[.6 .6 .6])
line([-3 2],[0 0],'color','k')
xlim([-3 2])
ylim([-1 4])
xlabel('Time (s)')
ylabel('Selectivity Contra-Ipsi (spk/s)')





%% population response correlation time course
PSTH_yes_tmp = PSTH_yes_cue_aligned;
PSTH_no_tmp = PSTH_no_cue_aligned;
for i_cell = 1:size(PSTH_yes_tmp,1)
    PSTH_yes_tmp(i_cell,:) = PSTH_yes_tmp(i_cell,:)-mean(PSTH_yes_tmp(i_cell,:));
    PSTH_no_tmp(i_cell,:) = PSTH_no_tmp(i_cell,:)-mean(PSTH_no_tmp(i_cell,:));

    PSTH_yes_tmp(i_cell,:) = PSTH_yes_tmp(i_cell,:)/var(PSTH_yes_tmp(i_cell,:));
    PSTH_no_tmp(i_cell,:) = PSTH_no_tmp(i_cell,:)/var(PSTH_no_tmp(i_cell,:));

end

i_sel = (cellType_all(:,1)==1 & sum(sig_selective,2)>0 & N_trials_all(:,1)>15 & N_trials_all(:,2)>15 & (spk_count_yes_all(:,4)>.1 & spk_count_no_all(:,4)>.1));

FR = [PSTH_yes_tmp(i_sel,:); PSTH_no_tmp(i_sel,:)];
FR_cue = FR(:,t==0);

Rpopulation_corr = [];
for i = 1:size(t,2)
    Rpopulation_corr(i,1) = corr(FR(:,i),FR_cue);
end

subplot(2,1,2);
plot(t,Rpopulation_corr);
line([0 0],[-.3 1],'color','k','linestyle',':')
line([0 0]-1.3,[-.3 1],'color','k','linestyle',':','color',[.6 .6 .6])
line([0 0]-2.6,[-.3 1],'color','k','linestyle',':','color',[.6 .6 .6])
line([-3 2],[0 0],'color','k')
xlim([-3 2])
ylim([-.3 1])
xlabel('Time (s)')
ylabel('Population response correlation')





%% cell type by selectivity

% all prymidal cells
cell_type_123 = zeros(size(sig_selective,1),1);
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & ~sig_selective(:,3) & cellType_all(:,1)==1)) = 1;
cell_type_123(find((sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==1)) = 2;
cell_type_123(find(~(sig_selective(:,1)|sig_selective(:,2)) & sig_selective(:,3) & cellType_all(:,1)==1)) = 3;
cell_type_tmp = [sum(cell_type_123==1) sum(cell_type_123==2) sum(cell_type_123==3)];
cell_type_btstrp = [];
for i_btstrp = 1:1000
    % sem over mice
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
    % sem over mice
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
subplot(1,2,1); hold on
h = bar(1, cell_type_tmp(1)/sum(cell_type_tmp));    set(h,'facecolor','w') 
h = bar(2, cell_type_tmp(2)/sum(cell_type_tmp));    set(h,'facecolor',[.8 .8 .8])
h = bar(3, cell_type_tmp(3)/sum(cell_type_tmp));    set(h,'facecolor','k') 
errorbar(1, cell_type_tmp(1)/sum(cell_type_tmp), std(cell_type_btstrp(:,1)),'k');
errorbar(2, cell_type_tmp(2)/sum(cell_type_tmp), std(cell_type_btstrp(:,2)),'k');
errorbar(3, cell_type_tmp(3)/sum(cell_type_tmp), std(cell_type_btstrp(:,3)),'k');
ylim([0 .6])
title('All prymidal neurons')
legend('preparatory','prep.+ peri','peri-movement')
ylabel('Fraction of neurons')

subplot(1,2,2); hold on
bar(5, preference(1)/sum(preference),'b');
bar(6, preference(2)/sum(preference),'r');
errorbar(5, preference(1)/sum(preference), std(preference_btstrp(:,1)),'k');
errorbar(6, preference(2)/sum(preference), std(preference_btstrp(:,2)),'k');
line([4 7],[.5 .5 ],'linestyle',':','color','k')
ylim([0 .8])
legend('contra-preferring','ipsi-preferring')





%% selectivity time course, binned
t =  -3:.1:1.8;
sig_all = [];
pref_all = [];

i_t = 0;
for i_time = -3:.1:1.8
    
    i_t = i_t+1;
    
    t_start = i_time-.1;
    t_end = i_time+.1;

    n_cell = 0;
    for i_cell = 1:size(spike_times_all,1)
        
        if N_trials_all(i_cell,1)>15 & N_trials_all(i_cell,2)>15  & cellType_all(i_cell,1)==1;
        
            n_cell = n_cell+1;
            
            spk_yes_tmp = [];
            for i_trial = 1:size(spike_times_all{i_cell,1});
                spk_t = spike_times_all{i_cell,1}{i_trial};
                spk_yes_tmp(i_trial,1) = sum(spk_t>t_start & spk_t<t_end);
            end
            spk_no_tmp = [];
            for i_trial = 1:size(spike_times_all{i_cell,2});
                spk_t = spike_times_all{i_cell,2}{i_trial};
                spk_no_tmp(i_trial,1) = sum(spk_t>t_start & spk_t<t_end);
            end
            sig_all(n_cell,i_t) = ttest2(spk_yes_tmp,spk_no_tmp);
            pref_all(n_cell,i_t) = mean(spk_yes_tmp) - mean(spk_no_tmp);
        end
    end
    
end


sig_all(isnan(sig_all))=0;

figure; hold on
bar(t, sum(sig_all & pref_all>0),'b')
bar(t, -sum(sig_all & pref_all<0),'r')
line([-1.3 -1.3],[-300 300],'color','k')
line([-1.3 -1.3]*2,[-300 300],'color','k')
line([0 0],[-300 300],'color','k')
xlim([-3 2])
ylabel('Number of neurons')








