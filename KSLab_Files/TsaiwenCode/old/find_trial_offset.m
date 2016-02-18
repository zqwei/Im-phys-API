function [offset,maxcorr]=find_trial_offset(trial_duration,nimage)
    nimage=nimage(1:end-1);
    nbehav_trials=length(trial_duration);
    nim_trials=length(nimage);
    for i=0:(nbehav_trials-nim_trials)
        trials=trial_duration((1:length(nimage))+i);
        cc=corrcoef(trials,nimage);
        coef(i+1)=cc(1,2);
    end
    figure;plot(coef);
    [maxcorr,offset]=max(coef);
    offset=offset-1;
    figure;plot(nimage,trial_duration((1:length(nimage))+offset),'.')
end