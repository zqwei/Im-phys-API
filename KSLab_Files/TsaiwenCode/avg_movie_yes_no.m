% bav_file=dir('../../../behav/alltrial.mat');


files=dir('Image*.tif');
lastfile=files(end).name;
ntrial=str2num(lastfile(end-6:end-4))-1;
prefix=lastfile(1:end-7);

nimage=zeros(1,ntrial);
for i=1:ntrial
    name=[prefix,num2str(i,'%03d'),'.tif'];
    a=dir(name);
    if isempty(a)
        nimage(i)=0;       
    else
        info=imfinfo(a.name);
        nimage(i)=length(info);
    end        
end
%%
load('../../../behav/alltrial.mat');
[offset,maxcorr]=find_trial_offset_handle_missing([alltrial.duration],nimage)

%%
bev_range=(1:length(nimage))+offset;
trial=alltrial(bev_range);
duration=[trial.duration];
% figure;plot(duration,nimage,'.');

ntrial=length(trial);
frameRate=29.68/2;   %15HZ
range=-45:45;
delay=25:45;
sample=7:24;
resp=46:60;

% frameRate=7.22;      % 7 Hz 4 planes
% range=-50:21;
% delay=5:50;
% sample=7:24;
% resp=46:60;

cueframe=floor(frameRate*[trial.cuetime]);
type=[trial.type];


correct=[trial.correct];
early=[trial.early];
idx1=(correct&(type=='r')&(~early))&(nimage>0);
idx2=(correct&(type=='l')&(~early))&(nimage>0);
disp(' ');
disp(['% correct:   ', num2str(sum(correct)/length(trial)*100)]);
disp(['# trials:    ', num2str(length(trial))]);
disp(['# correct:   ',num2str(sum(correct))]);
disp(['# early:     ',num2str(sum(early))]);
disp(['useful R:    ', num2str(sum(idx1))]);
disp(['useful L:    ', num2str(sum(idx2))]);

% idx1(1:80)=0;
% idx2(1:80)=0;
%%
rtrials=find(idx1);
avg_r=zeros(512,512,length(range));
for i=1:length(rtrials)
    name=[prefix,num2str(rtrials(i),'%03d'),'.tif'];
    r=cueframe(rtrials(i))+range;
    [im hdr] = load_image(name,[min(r),max(r)]);
%     avg_r=avg_r+im(:,:,cueframe(rtrials(i))+range);
    avg_r=avg_r+im;
    i
end
avg_r=avg_r/length(rtrials);

%%
ltrials=find(idx2);
avg_l=zeros(512,512,length(range));
for i=1:length(ltrials)
    name=[prefix,num2str(ltrials(i),'%03d'),'.tif'];
    r=cueframe(ltrials(i))+range;
    [im hdr] = load_image(name,[min(r),max(r)]);
%     avg_r=avg_r+im(:,:,cueframe(rtrials(i))+range);
    avg_l=avg_l+im;
    i
end
avg_l=avg_l/length(ltrials);


%% merge movie
bin_r=avg_r(1:2:end,1:2:end,:)+avg_r(1:2:end,2:2:end,:)+avg_r(2:2:end,1:2:end,:)+avg_r(2:2:end,2:2:end,:);
bin_l=avg_l(1:2:end,1:2:end,:)+avg_l(1:2:end,2:2:end,:)+avg_l(2:2:end,1:2:end,:)+avg_l(2:2:end,2:2:end,:);


base=1:5;


f0_r=mean(bin_r(:,:,base),3)+1;
f0_l=mean(bin_l(:,:,base),3)+1;

f0_r=repmat(f0_r,[1,1,length(range)]);
f0_l=repmat(f0_l,[1,1,length(range)]);
dff_r=(bin_r-f0_r)./f0_r;
dff_l=(bin_l-f0_l)./f0_l;
ss=size(dff_r);
merge=zeros([ss,3]);
merge(:,:,:,1)=dff_r;
merge(:,:,:,2)=dff_l;
merge(merge>2)=2;
merge(merge<-2)=-2;
merge=merge/2;

% merge(1:30,1:30,sample,1)=1;
% merge(1:30,1:30,sample,2)=0;
% merge(1:30,1:30,sample,3)=0;
% 
% merge(1:30,1:30,delay,1)=0;
% merge(1:30,1:30,delay,2)=1;
% merge(1:30,1:30,delay,3)=1;


matVis(merge);
save('avg_summary','avg_r','avg_l','merge');
% %%
% % base=1:5;
% % resp=14:21;
% 
% base=1:5;
% % r=sample;
% r=resp;
% f0_r=mean(avg_r(:,:,base),3);
% f0_l=mean(avg_l(:,:,base),3);
% dff_r=(mean(avg_r(:,:,r),3)-f0_r)./f0_r;
% dff_l=(mean(avg_l(:,:,r),3)-f0_l)./f0_l;
% 
% dff=cat(2,dff_r,dff_l);
% % f0=mean(c(:,:,base),3);
% % df=mean(c(:,:,resp),3)-f0;
% 
% 
% ss=size(dff_r);
% merge=zeros([ss(1:2),3]);
% merge(:,:,1)=dff_r/2;
% merge(:,:,2)=dff_l/2;
% merge(merge>1)=1;
% merge(merge<0)=0;
% figure;image(merge);

%%
nimage=size(merge,3);
hf=figure('position',[100,100,256,256]);
h=image(zeros(256,256,3));
set(gca,'position',[0,0,1,1],'visible','off');
axis image;

htext=text(10,10,'baseline','color','w','FontSize',20);
% axis image;
% drawnow;
ax=gca;
% set(ax,'visible','off');

aviobj = avifile('example.avi','compression','None','fps',7);
pause(2);

for i=1:nimage
%     test=im(:,:,i);
%     RGB=repmat(test,[1,1,3]);
%     RGB=(RGB)/(75);
    RGB=squeeze(merge(:,:,i,:));
    RGB(RGB>1)=1;
    RGB(RGB<0)=0;
    set(h,'CData',RGB);
%    
    F=getframe(ax);
    aviobj=addframe(aviobj,F);
    if ismember(i,sample)
         set(htext,'String','sample');
    elseif ismember(i,delay)
         set(htext,'String','delay');
    elseif i>45
         set(htext,'String','response');
    end
        
end
aviobj=close(aviobj);

