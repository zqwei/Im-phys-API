% file=dir('data*.mat');
file=dir('an*.mat');
load(file(1).name);

% #43 Correct right
% #44 Incorrect right
% #45 No response right
% #51 Correct left
% #52 Incorrect left
% #53 No response left
% #56/57 alarm

event=saved_history.RewardsSection_LastTrialEvents;
duration=[];
trialtype=(saved.yes_no_multi_pole_twcobj_trial_type);  %'r'=114; 'l'=108
alltrial=[];
for i=2:length(event)
    trial=[];
    trig_idx=find(event{i}(:,1)==40);
    start_time=event{i}(trig_idx(1),3);
    if i<length(event)
        end_time=event{i+1}(trig_idx(1),3);        
    else
        end_time=event{i}(end,3);
    end
    trial.duration=end_time-start_time;
        
%     start_time=event{i}(1,3);           
    trial.type=char(trialtype(i)); 
    trial.correct=(sum(event{i}(:,1)==43)+sum(event{i}(:,1)==51))>0;
    trial.noresponse=(sum(event{i}(:,1)==45)+sum(event{i}(:,1)==53))>0;
    trial.early=(sum(event{i}(:,1)==56)+sum(event{i}(:,1)==57))>0;
    

%     trial.duration=event{i}(end,3)-start_time;
    trial.leftlicktime=event{i}((event{i}(:,2)==1),3)-start_time;
    trial.rightlicktime=event{i}((event{i}(:,2)==3),3)-start_time;
    idx=find(event{i}(:,1)==55);  %%cue
    trial.cuetime=event{i}(idx(1),3)-start_time;
    alltrial=[alltrial,trial];    
end

correct=[alltrial.correct];
type=[alltrial.type];
early=[alltrial.early];

disp(['Performance: ',num2str(sum(correct)/length(correct))]);
disp(['Early resp:  ',num2str(sum(early)/length(early))]);
save('alltrial','alltrial');
%%

subtrial=alltrial(correct&(type=='l')&(~early));

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
figure;hold on;

for i=1:length(alltrial)
    if alltrial(i).correct
        if ~(alltrial(i).early)
            color='g.';
        else
            color='b.';
        end
    else
        color='r.';
    end
    plot(i,alltrial(i).type=='l',color);
end
ylim([-1,2]);