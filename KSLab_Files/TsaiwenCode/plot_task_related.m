function [dff_avg,pval,select,ROI_list]=plot_task_related(file,plotflag)
if ~exist('file')
    file=uigetfile('*para.mat');
    plotflag=1;
end
load(file);
% load('an022_2013_08_23_580_para.mat')
%%
% idx=strcmp({ROI_list.type},'p');
% ROI_list=ROI_list(idx);
%%
nROI=length(ROI_list);
range=-45:45;
frameRate=29.68/2;
ntgroup=18;
nbin=5;

t=range/frameRate;
dff_avg=zeros(ntgroup,nROI,2);
pval=zeros(nROI,2);
stat={};
signific=zeros(ntgroup,nROI,2);
select=zeros(ntgroup,nROI);
sindex=zeros(ntgroup,nROI);
for i=1:length(ROI_list)
    [falign,idx1,idx2]=plot_align_response(ROI_list(i),trial(1:end-1),nimage(1:end-1),range,0);
    falign=falign(1:(nbin*ntgroup),:);
    falign=squeeze(mean(reshape(falign,nbin,ntgroup,[])));
    [pval(i,1),table,stat{i,1}]= kruskalwallis(falign(:,idx1)',[],'off');    
    [pval(i,2),table,stat{i,2}]= kruskalwallis(falign(:,idx2)',[],'off');      

    dff_avg(:,i,1)=rawf2df_f(median(falign(:,idx1),2),1:2);    
    dff_avg(:,i,2)=rawf2df_f(median(falign(:,idx2),2),1:2);
    dff_align=rawf2df_f(falign,1);
    for k=1:ntgroup
        p=ranksum(dff_align(k,idx1),dff_align(k,idx2));
        select(k,i)=p<0.01;        
        sindex(k,i)=(dff_avg(k,i,1)-dff_avg(k,i,2))/(abs(dff_avg(k,i,1))+abs(dff_avg(k,i,2)));
    end
        
end

%%
if plotflag
idx=find(sum((pval<0.001),2)>0);
figure;subplot(1,2,1);imagesc(dff_avg(:,idx,1)');set(gca,'clim',[0,1]);
subplot(1,2,2);imagesc(dff_avg(:,idx,2)');set(gca,'clim',[0,1]);

figure;imagesc(select(:,idx)'.*sindex(:,idx)');set(gca,'clim',[-1,1]);
end
select=select.*sindex;