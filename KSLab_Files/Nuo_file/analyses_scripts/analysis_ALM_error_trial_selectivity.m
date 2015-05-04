%% supplemental figure 2


clear all
close all

% check to see if the data has been compiled
if ~exist('ALM_compiled_all_data.mat')
    
    % if not, compile all data first
    analysis_compile_all_data
    
end


load ALM_compiled_all_data

sig_selective(isnan(sig_selective))=0;







%% example neurons, PSTH

example_cells = {
    'ANM218693', 3, 8;...
    'ANM219247', 2, 4;...
    'ANM218693', 1, 4;...
    'ANM221977', 2, 8;...
    'ANM218457', 4, 9;...
    };

figure
n_plot = 0;
for i_ex_cell = 1:size(example_cells,1)
    i_mice = find(strcmp(animal_ID,example_cells{i_ex_cell,1}));
    i_session = example_cells{i_ex_cell,2};
    i_unit = example_cells{i_ex_cell,3};
    n_plot = n_plot+1;
    
    i_cell = find(Mice_all(:,1)==i_mice & Mice_all(:,2)==i_session & Mice_all(:,3)==i_unit);

    

    subplot(2,5,n_plot); hold on
    plot(t,PSTH_yes_cue_aligned(i_cell,:),'color','b');
    plot(t,PSTH_no_cue_aligned(i_cell,:),'color','r');
    
    ymax = max([PSTH_yes_cue_aligned(i_cell,:) PSTH_no_cue_aligned(i_cell,:) PSTH_yes_error_cue_aligned(i_cell,:) PSTH_no_error_cue_aligned(i_cell,:)]);
    line([0 0],[0 ymax*1.2],'color','k')
    line([-1.3 -1.3],[0 ymax*1.2],'color','k')
    line([-2.6 -2.6],[0 ymax*1.2],'color','k')
    ylim([0 ymax*1.25])
    xlim([-3.1 2]);
    
    
    subplot(2,5,n_plot+5); hold on
    plot(t,PSTH_yes_error_cue_aligned(i_cell,:),'color',[.7 .7 1]);
    plot(t,PSTH_no_error_cue_aligned(i_cell,:),'color',[1 .7 .7]);
    
    ymax = max([PSTH_yes_cue_aligned(i_cell,:) PSTH_no_cue_aligned(i_cell,:) PSTH_yes_error_cue_aligned(i_cell,:) PSTH_no_error_cue_aligned(i_cell,:)]);
    line([0 0],[0 ymax*1.2],'color','k')
    line([-1.3 -1.3],[0 ymax*1.2],'color','k')
    line([-2.6 -2.6],[0 ymax*1.2],'color','k')
    ylim([0 ymax*1.25])
    xlim([-3.1 2]);
    
    
    if n_plot==1
        subplot(2,5,n_plot); hold on
        xlabel('Time (s)')
        ylabel('Spk/s')
        title('correct trials')
        legend('lick right','lick left')

    
        subplot(2,5,n_plot+5); hold on
        title('error trials')
    
    end
    

end











%% correct vs. error trials, scatter
i_selective = (sig_selective(:,1)|sig_selective(:,2)|sig_selective(:,3)) & cellType_all(:,1)==1;
i_n_trial = N_trials_all(:,1)>10 & N_trials_all(:,2)>10 & N_trials_all(:,3)>5 & N_trials_all(:,4)>5;

trial_selectivity = spk_count_yes_all-spk_count_no_all;
error_trial_selectivity = spk_count_yes_error_all-spk_count_no_error_all;
% % i_selective = i_selective & ((abs(trial_selectivity(:,1)-error_trial_selectivity(:,1))<4) & abs(trial_selectivity(:,1))>2);% | ((abs(trial_selectivity(:,1)-error_trial_selectivity(:,1))<1) & abs(trial_selectivity(:,1)<2));

trial_selectivity(trial_selectivity>15) = 15;
trial_selectivity(trial_selectivity<-15) = -15;
error_trial_selectivity(error_trial_selectivity>15) = 15;
error_trial_selectivity(error_trial_selectivity<-15) = -15;


figure;
subplot(1,3,1); hold on
plot(trial_selectivity(i_n_trial,1),error_trial_selectivity(i_n_trial,1),'.','color',[.7 .7 .7]);%'ok');
plot(trial_selectivity(i_n_trial & i_selective,1),error_trial_selectivity(i_n_trial & i_selective,1),'.k');%,'ok','markerfacecolor','k');
line([-15 15],[-15 15],'color','k');
line([-15 15],[15 -15],'color','k');
line([-15 15],[0 0],'color','k');
line([0 0],[-15 15],'color','k');
xlim([-15 15])
ylim([-15 15])
disp('***** sample *****')
[r p] = corr(trial_selectivity(i_n_trial & i_selective,1),error_trial_selectivity(i_n_trial & i_selective,1))

subplot(1,3,2); hold on
plot(trial_selectivity(i_n_trial,2),error_trial_selectivity(i_n_trial,2),'.','color',[.7 .7 .7]);%'ok');
plot(trial_selectivity(i_n_trial & i_selective,2),error_trial_selectivity(i_n_trial & i_selective,2),'.k');%'ok','markerfacecolor','k');
line([-15 15],[-15 15],'color','k');
line([-15 15],[15 -15],'color','k');
line([-15 15],[0 0],'color','k');
line([0 0],[-15 15],'color','k');
xlim([-15 15])
ylim([-15 15])
disp('***** delay *****')
[r p] = corr(trial_selectivity(i_n_trial & i_selective,2),error_trial_selectivity(i_n_trial & i_selective,2))


subplot(1,3,3); hold on
plot(trial_selectivity(i_n_trial,3),error_trial_selectivity(i_n_trial,3),'.','color',[.7 .7 .7]);%'ok');
plot(trial_selectivity(i_n_trial & i_selective,3),error_trial_selectivity(i_n_trial & i_selective,3),'.k');%'ok','markerfacecolor','k');
line([-15 15],[-15 15],'color','k');
line([-15 15],[15 -15],'color','k');
line([-15 15],[0 0],'color','k');
line([0 0],[-15 15],'color','k');
xlim([-15 15])
ylim([-15 15])
disp('***** response *****')
[r p] = corr(trial_selectivity(i_n_trial & i_selective,3),error_trial_selectivity(i_n_trial & i_selective,3))






