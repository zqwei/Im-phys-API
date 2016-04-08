file=dir('*.mat');
load(file(1).name);

% #54 Cue

event=saved_history.RewardsSection_LastTrialEvents;
duration=[];
trialtype=(saved.yes_no_multi_pole_twcobj_trial_type);  %'r'=114; 'l'=108
alltrial=[];
for i=2:length(event)
    trial=[];
    trig_idx=find(event{i}(:,1)==40);
    start_time=event{i}(trig_idx(1),3);

    trial.type=char(trialtype(i));   

    trial.duration=event{i}(end,3)-start_time;
    trial.leftlicktime=event{i}((event{i}(:,2)==1),3)-start_time;
    trial.rightlicktime=event{i}((event{i}(:,2)==3),3)-start_time;
    idx=find(event{i}(:,1)==54);  %%cue
    trial.cuetime=event{i}(idx(1),3)-start_time;
    alltrial=[alltrial,trial];    
end


%%
type=[alltrial.type];
subtrial=alltrial((type=='r'));

figure;hold on;
for i=1:length(subtrial)
    if ~isempty(subtrial(i).leftlicktime)
        plot(subtrial(i).leftlicktime-subtrial(i).cuetime,i,'.r');
    end
    if ~isempty(subtrial(i).rightlicktime)
        plot(subtrial(i).rightlicktime-subtrial(i).cuetime,i+0.2,'.k');    
    end
end
xlim([-4,8]);

%%
% figure;hold on;
% 
% for i=1:length(alltrial)
%     if alltrial(i).correct
%         if ~(alltrial(i).early)
%             color='g.';
%         else
%             color='b.';
%         end
%     else
%         color='r.';
%     end
%     plot(i,alltrial(i).type=='l',color);
% end
% ylim([-1,2]);