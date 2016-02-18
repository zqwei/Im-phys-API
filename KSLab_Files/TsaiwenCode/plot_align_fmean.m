function [falign,idx1,idx2]=plot_align_fmean(fmean,trial,nimage,plot_flag,overlay_red,overlay_farred)
frameRate=29.68/2;
xrange=-45:45;


% frameRate=29.68/4;
% xrange=-60:20;
% polein=-7.2*frameRate;
% poleout=-6*frameRate;
% cue=0;


% bev_range=(1:92)+110;
% im_range=1:92;
% figure;plot([alltrial(bev_range).duration],nimage(im_range),'.');

% fmean=fmean_neuropil;
valid=nimage>0;

%%
% trial=alltrial(bev_range);
ntrial=length(trial);
% fmean=ROI_list(1).fmean;

cueframe=floor(frameRate*[trial.cuetime]);
offset=[0,cumsum(nimage(1:(end-1)))];


cueframe=cueframe+offset;

temp=fmean;temp(temp==0)=[];
f0=mode(temp);
dff=(fmean-f0)/f0;
falign=zeros(length(xrange),ntrial);
for i=1:length(cueframe)
    falign(:,i)=fmean(xrange+cueframe(i));
%     dffalign(:,i)=dff(xrange+cueframe(i));
    dffalign(:,i)=rawf2df_f(falign(:,i),1:5);
end

%%

correct=[trial.correct];
type=[trial.type];
early=[trial.early];
idx1=(correct&(type=='r')&(~early))&(valid);
idx2=(correct&(type=='l')&(~early))&(valid);
t=xrange/frameRate;

%%
if plot_flag
yrange=[-0.2,0.5];

figure('position',[10,500,600,250]);
hold on;
plot(t,dffalign(:,idx1),'color',[1,0.8,0.8]);
plot(t,dffalign(:,idx2),'color',[0.8,1,0.8]);
h1=plot(t,mean(dffalign(:,idx1),2),'r','linewidth',3);
h2=plot(t,mean(dffalign(:,idx2),2),'g','linewidth',3);
plot([-2.6,-2.6],[-1,4],'b:','linewidth',2);
plot([-1.4,-1.4],[-1,4],'b:','linewidth',2);
plot([0,0],[-1,4],'b:','linewidth',2);
ylim(yrange)
xlim([-8,4.5]);
box off;
legend([h1,h2],'correct right','correct left');

% figure;
% subplot(1,2,1);
% % image(overlay_red);
% axis image;
% hold on;
% plot(ROI.centerPos(2),ROI.centerPos(1),'+w');
% subplot(1,2,2);
% image(overlay_farred);
% axis image;hold on;
% plot(ROI.centerPos(2),ROI.centerPos(1),'+w');
figure('position',[10,500,600,250]);

subplot(1,2,1);
title('correct right');
plot(t,dffalign(:,idx1),'color',[1,1,1]*0.8);
hold on;
plot(t,mean(dffalign(:,idx1),2),'r','linewidth',3);
plot([-2.6,-2.6],[-1,4],'b:','linewidth',2);
plot([-1.4,-1.4],[-1,4],'b:','linewidth',2);
plot([0,0],[-1,4],'b:','linewidth',2);

ylim(yrange)
xlim([-4,4.5]);
box off;
subplot(1,2,2);
title('correct left');
hold on;
plot(t,dffalign(:,idx2),'color',[1,1,1]*0.8);
plot(t,mean(dffalign(:,idx2),2),'g','linewidth',3);


plot([-2.6,-2.6],[-1,4],'b:','linewidth',2);
plot([-1.4,-1.4],[-1,4],'b:','linewidth',2);
plot([0,0],[-1,4],'b:','linewidth',2);
ylim(yrange)
xlim([-4,4.5]);
box off;

figure;
subplot(1,2,1);hold on;
imagesc(t,[1,sum(idx1)],dffalign(:,idx1)')
set(gca,'clim',[yrange]);
plot([-2.6,-2.6],[0,sum(idx1)],'w:','linewidth',2);
plot([-1.4,-1.4],[0,sum(idx1)],'w:','linewidth',2);
plot([0,0],[0,sum(idx1)],'w:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,sum(idx1)]);
title('correct right');
box off;


subplot(1,2,2);hold on;
imagesc(t,[1,sum(idx2)],dffalign(:,idx2)')
set(gca,'clim',[yrange]);
plot([-2.6,-2.6],[0,sum(idx2)],'w:','linewidth',2);
plot([-1.4,-1.4],[0,sum(idx2)],'w:','linewidth',2);
plot([0,0],[0,sum(idx2)],'w:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,sum(idx2)]);
title('correct left');
box off;


% figure;imagesc(mask);

end