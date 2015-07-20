load('an026_2013_11_12_545_para.mat','trial','nimage','ROI_list');
for i=1:length(trial)
    all_lick=[trial(i).leftlicktime;trial(i).rightlicktime];
    t=[sort(all_lick);trial(i).duration];
    isolate=[];
    for k=2:(length(t)-1)
    	if ((t(k)-t(k-1))>1)&&((t(k+1)-t(k))>1)
            isolate=[isolate,t(k)];
        end
    end
    trial(i).isolate=isolate;
    trial(i).exist_iso=~isempty(isolate);
end