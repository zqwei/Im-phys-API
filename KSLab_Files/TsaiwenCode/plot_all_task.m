%% pile data
files=uigetfile('*para.mat','multiselect','on');
dff_avg=[];
pval=[];
select=[];
ROI_list=[];

for i=1:length(files)
%     a=plot_population(files{i},0);
    
    
    [dff,p,s,roi]=plot_task_related(files{i},0);
    dff_avg=cat(2,dff_avg,dff);   
    pval=cat(1,pval,p);
    select=cat(2,select,s);
    ROI_list=[ROI_list,roi];
    i
end

save('pile_result','ROI_list','pval','select','dff_avg');
%%
depth=[];
for i=1:length(files)
    s=load(files{i},'depth','ROI_list');
    depth=[depth;ones(length(s.ROI_list),1)*s.depth];
end
save('pile_result','depth','-append');
%%
load pile_result;
fm=zeros(size(ROI_list));
fm_neuropil=zeros(size(ROI_list));
for i=1:length(ROI_list)
    fm(i)=mean(ROI_list(i).fmean);
    fm_neuropil(i)=mean(ROI_list(i).fmean_neuropil);
end
inc=(fm./fm_neuropil>1.03);
%%
% late
% sel_right=(abs(sum(select(2:9,:)))<2)'&(sum(select(10:14,:))>2)';
% sel_left=(abs(sum(select(2:9,:)))<2)'&(sum(select(10:14,:))<-2)';
% 
% % early 
% sel_right=(sum(select(2:9,:))>2)';
% sel_left=(sum(select(2:9,:))<-2)';

% % all period
sel_right=(sum(select(2:18,:))>3)';
sel_left=(sum(select(2:18,:))<-3)';

pt=(strcmp({ROI_list.type},'p'))';
it=(strcmp({ROI_list.type},'i'))';
task=sum((pval<0.001),2)>0;

%all stat:
disp('all cell statistics ');
disp(['  # all cell     : ',num2str(sum(inc))]);
disp(['  # task related : ',num2str(sum(inc'&task))]);
disp(['  # contra pref  : ',num2str(sum(inc'&task&sel_right))]);
disp(['  # ipsi pref    : ',num2str(sum(inc'&task&sel_left))]);
disp(' ');

%%
%deep stat:
disp('Deep cell statistics (z>450um)');
disp(['  # all cell     : ',num2str(sum((depth>450)&inc'))]);
disp(['  # task related : ',num2str(sum((depth>450)&inc'&task))]);
disp(['  # contra pref  : ',num2str(sum((depth>450)&inc'&task&sel_right))]);
disp(['  # ipsi pref    : ',num2str(sum((depth>450)&inc'&task&sel_left))]);
disp(' ');

disp('Deep PT cell statistics (z>450um)');
disp(['  # all cell     : ',num2str(sum((depth>450)&inc'&pt))]);
disp(['  # task related : ',num2str(sum((depth>450)&inc'&task&pt))]);
disp(['  # contra pref  : ',num2str(sum((depth>450)&inc'&task&sel_right&pt))]);
disp(['  # ipsi pref    : ',num2str(sum((depth>450)&inc'&task&sel_left&pt))]);
disp(' ');

disp('Deep IT cell statistics (z>450um)');
disp(['  # all cell     : ',num2str(sum((depth>450)&inc'&it))]);
disp(['  # task related : ',num2str(sum((depth>450)&inc'&task&it))]);
disp(['  # contra pref  : ',num2str(sum((depth>450)&inc'&task&sel_right&it))]);
disp(['  # ipsi pref    : ',num2str(sum((depth>450)&inc'&task&sel_left&it))]);
disp(' ');

disp('Superficial cell statistics (z<350um)');
disp(['  # all cell     : ',num2str(sum((depth<350)&inc'))]);
disp(['  # task related : ',num2str(sum((depth<350)&inc'&task))]);
disp(['  # contra pref  : ',num2str(sum((depth<350)&inc'&task&sel_right))]);
disp(['  # ipsi pref    : ',num2str(sum((depth<350)&inc'&task&sel_left))]);
disp(' ');

disp('Superficial IT cell statistics (z<350um)');
disp(['  # all cell     : ',num2str(sum((depth<350)&inc'&it))]);
disp(['  # task related : ',num2str(sum((depth<350)&inc'&task&it))]);
disp(['  # contra pref  : ',num2str(sum((depth<350)&inc'&task&sel_right&it))]);
disp(['  # ipsi pref    : ',num2str(sum((depth<350)&inc'&task&sel_left&it))]);
disp(' ');

%%
type={ROI_list.type};
% type_idx=strcmp(type,'p');
% type_idx=~(strcmp(type,'p')|strcmp(type,'i'));
type_idx=1;
task=sum((pval<0.001),2)>0;
% sel=task&type_idx'&(depth>450)&inc';
sel=task&inc';
dff=dff_avg(:,sel,:);
select=select(:,sel);
%% calculate first frame index
si=sum(select);
first_ind=zeros(length(si),1);
for i=1:length(first_ind)
    if si(i)>0
        temp=find(select(:,i)>0);
        first_ind(i)=temp(1);
    elseif si(i)<0
        temp=find(select(:,i)<0);
        first_ind(i)=temp(1);
    end
end
%%

frameRate=29.68/2;
t=(-45:44)/frameRate;
t=mean(reshape(t,5,[]));
resp=3:12;

[h,idx1]=sort(mean(dff(resp,:,1),1),'descend');
[h,idx2]=sort(mean(dff(resp,:,2),1),'descend');

plot_multi_roi(t,dff(:,idx1,1),dff(:,idx2,2));

%% right preferring 
h=sum(select);
[x,ind]=sort(first_ind);
h=h(ind);
temp=dff(:,ind,:);
idx=h>3;

plot_multi_roi(t,temp(:,idx,1),temp(:,idx,2));
subplot(1,4,4);
y=hist(x(idx),1:18);
bar(t,y,'facecolor',[1,1,1]*0.7,'edgecolor','none');hold on;
plot([-2.6,-2.6],[-0.2,max(y)*1.3],'k:','linewidth',2);
plot([-1.4,-1.4],[-0.2,max(y)*1.3],'k:','linewidth',2);
plot([0,0],[-0.2,max(y)*1.3],'k:','linewidth',2);box off;
xlim([min(t),max(t)]);
box off;
sum(y(1:9))
sum(y(10:end))
disp(['# right preferring: ',num2str(sum(idx))]);
%% left preferring 
h=sum(select);
[x,ind]=sort(first_ind);
h=h(ind);
temp=dff(:,ind,:);
idx=h<-3;
plot_multi_roi(t,temp(:,idx,1),temp(:,idx,2));
subplot(1,4,4);
y=hist(x(idx),1:18);
bar(t,y,'facecolor',[1,1,1]*0.7,'edgecolor','none');hold on;
plot([-2.6,-2.6],[-0.2,max(y)*1.3],'k:','linewidth',2);
plot([-1.4,-1.4],[-0.2,max(y)*1.3],'k:','linewidth',2);
plot([0,0],[-0.2,max(y)*1.3],'k:','linewidth',2);box off;
xlim([min(t),max(t)]);


box off;
sum(y(1:9))
sum(y(10:end))
disp(['# left preferring: ',num2str(sum(idx))]);
%%
[h,idx]=sort(sum(select),'descend');
% figure;imagesc(select(:,idx)');
x=1:length(idx);
figure;
subplot(1,3,1:2);hold on;
imagesc(t,[0,length(idx)],-select(:,idx)');set(gca,'clim',[-1,1]);
plot([-2.6,-2.6],[0,length(idx)],'k:','linewidth',2);
plot([-1.4,-1.4],[0,length(idx)],'k:','linewidth',2);
plot([0,0],[0,length(idx)],'k:','linewidth',2);
xlim([min(t),max(t)]);ylim([0,length(idx)]);
set(gca,'YDir','reverse');
subplot(1,3,3);hold on;

plot(h,x,'k','linewidth',2);set(gca,'YDir','reverse');
plot(-3*[1,1],[0,length(idx)],':');
plot(3*[1,1],[0,length(idx)],':')
ylim([0,length(idx)]);