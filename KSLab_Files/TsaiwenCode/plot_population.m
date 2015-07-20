% file=uigetfile('*para.mat');
% plotflag=1;

function [dff_avg,idx1,idx2]=plot_population(file,plotflag)
if ~exist('file')
    file=uigetfile('*para.mat');
    plotflag=1;
end
load(file,'ROI_list','trial','nimage');
% %%
idx=strcmp({ROI_list.type},'p');
% ROI_list=ROI_list(idx);
%%
frameRate=29.68/2;
% resp=30:60;

% range=-70:45;
% polein=-4.2;
% poleout=-3;

range=-45:45;
polein=-2.6;
poleout=-1.4;




%% delay 1.4s
% frameRate=29.68/4;
% range=-15:45;

% %% delay 3s
% frameRate=29.68/4;
% range=-35:45;
% polein=-4.2;
% poleout=-3;
% 
%% delay 6s
% frameRate=29.68/4;
% range=-60:45;
% polein=-7.2;
% poleout=-6;
%%

t=range/frameRate;
dff_avg=zeros(length(range),length(ROI_list),2);

peak_loc=zeros(length(ROI_list),1);
peak_size=zeros(length(ROI_list),1);

for i=1:length(ROI_list)
    [falign,idx1,idx2]=plot_align_response(ROI_list(i),trial(1:end-1),nimage(1:end-1),range,0);
    dff_avg(:,i,1)=median(rawf2df_f((falign(:,idx1)),1:5),2);    
    dff_avg(:,i,2)=median(rawf2df_f((falign(:,idx2)),1:5),2);
    
    [M,idx]=max(mean(dff_avg(:,i,1),3));
    peak_loc(i)=idx;
    peak_size(i)=M;
    
    [M,idx]=max(mean(dff_avg(:,i,1),3));
    peak_loc(i)=idx;
end


% [s,idx1]=sort(mean(dff_avg(resp,:,1),1));
% [s,idx2]=sort(mean(dff_avg(resp,:,2),1));
% ROI_list=ROI_list(idx1);
idx=(1:length(peak_size))';

%% bin
% ss=size(dff_avg);
% dff_avg=dff_avg(1:90,:,:);
% dff_avg=squeeze(mean(reshape(dff_avg,3,30,[],2)));
% t=t(1:3:end);


%%
if plotflag==1

    
figure;
subplot(1,2,1); hold on;
imagesc(t,[1,size(idx,1)],dff_avg(:,idx,1)');
set(gca,'clim',[0,1])

plot(polein*[1,1],[0,size(idx,1)],'w:','linewidth',2);
plot(poleout*[1,1],[0,size(idx,1)],'w:','linewidth',2);
plot([0,0],[0,size(idx,1)],'w:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,size(idx,1)]);

subplot(1,2,2); hold on;
imagesc(t,[1,size(idx1,1)],dff_avg(:,idx,2)');
set(gca,'clim',[0,1])

plot(polein*[1,1],[0,size(idx,1)],'w:','linewidth',2);
plot(poleout*[1,1],[0,size(idx,1)],'w:','linewidth',2);
plot([0,0],[0,size(idx,1)],'w:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,size(idx,1)]);

figure;plot(mean(dff_avg(:,idx,1),2),'b');hold on;plot(mean(dff_avg(:,idx,2),2),'r')
end
