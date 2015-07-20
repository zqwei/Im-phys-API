files=uigetfile('*para.mat','multiselect','on');

dff_avg=[];
for i=1:length(files)
    a=plot_population(files{i},0);
    dff_avg=cat(2,dff_avg,a);   
end
    
%%
range=-45:45;
frameRate=29.68/2;
t=range/frameRate;
% idx=(1:size(dff_avg,2))';
resp=30:60;
% [s,idx]=sort(max(dff_avg(resp,:,1),[],1));
[s,idx]=sort(mean(dff_avg(resp,:,1),1));
[s,idx1]=sort(mean(dff_avg(resp,:,2),1));
figure;
subplot(1,2,1); hold on;
imagesc(t,[1,size(idx,1)],dff_avg(:,idx,1)');
set(gca,'clim',[-0.2,1])

plot([-2.6,-2.6],[0,length(idx)],'w:','linewidth',2);
plot([-1.4,-1.4],[0,length(idx)],'w:','linewidth',2);
plot([0,0],[0,length(idx)],'w:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,length(idx)]);

subplot(1,2,2); hold on;
imagesc(t,[1,size(idx,1)],dff_avg(:,idx1,2)');
set(gca,'clim',[-0.2,1])

plot([-2.6,-2.6],[0,length(idx)],'w:','linewidth',2);
plot([-1.4,-1.4],[0,length(idx)],'w:','linewidth',2);
plot([0,0],[0,length(idx)],'w:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,length(idx)]);

