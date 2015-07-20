function [offset,maxcorr]=find_trial_offset_handle_missing(trial_duration,nimage)

valid_trial=false(size(nimage));
for i=1:(length(nimage)-1)
    if nimage(i)>0 && nimage(i+1)>0
        valid_trial(i)=1;
    end
end
nimage=nimage(1:end-1);
valid_trial=valid_trial(1:end-1);

nbehav_trials=length(trial_duration);

nim_trials=length(nimage);
for i=0:(nbehav_trials-nim_trials)
    trials=trial_duration((1:length(nimage))+i);
    trials=trials(valid_trial);
    nn=nimage(valid_trial)';
    idx=(trials<15);
    
    cc=corrcoef(trials(idx),nn(idx));
    coef(i+1)=cc(1,2);
end
figure;plot(coef);
[maxcorr,offset]=max(coef);
offset=offset-1;
trials=trial_duration((1:length(nimage))+offset);
trials=trials(valid_trial);
figure;plot(nn,trials,'.')