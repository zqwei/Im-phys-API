files=dir('*main*.tif');
files=files(2:end-1);
nimage=[];
for i=1:length(files)
    info = imfinfo(files(i).name);
    nimage(i)=length(info);
end
    
%%
offset=find_trial_offset([alltrial.duration],nimage);

%%

bev_range=(1:length(nimage))+offset;
trial=alltrial(bev_range);
duration=[trial.duration];
figure;plot(duration,nimage,'.');

ntrial=length(trial);
% fmean=ROI_list(1).fmean;
frameRate=29.68;
hdr(1).frameRate=hdr(1).frameRate/2;
cueframe=floor(hdr(1).frameRate*[trial.cuetime]);
type=[trial.type];


correct=[trial.correct];
early=[trial.early];
idx1=(correct&(type=='r')&(~early));
idx2=(correct&(type=='l')&(~early));

% idx1=(type=='r');
% idx2=(type=='l');

%%
range=-45:45;
% range=-29:45;


%%

rtrials=find(idx1);
avg_r=zeros(512,512,length(range));
for i=1:length(rtrials)
    [im hdr] = load_image(files(rtrials(i)).name);
    avg_r=avg_r+im(:,:,cueframe(rtrials(i))+range);
    i
end
avg_r=avg_r/length(rtrials);
%%
ltrials=find(idx2);
avg_l=zeros(512,512,length(range));
for i=1:length(ltrials)
    [im hdr] = load_image(files(ltrials(i)).name);
    avg_l=avg_l+im(:,:,cueframe(ltrials(i))+range);
    i
end
avg_l=avg_l/length(ltrials);

%%
base=1:15;
% resp=15:30;
% resp=31:45;
% resp=45:60;
resp=61:75;
f0_r=mean(avg_r(:,:,base),3);
f0_l=mean(avg_l(:,:,base),3);
dff_r=(mean(avg_r(:,:,resp),3)-f0_r)./f0_r;
dff_l=(mean(avg_l(:,:,resp),3)-f0_l)./f0_l;

dff=cat(2,dff_r,dff_l);
% f0=mean(c(:,:,base),3);
% df=mean(c(:,:,resp),3)-f0;
figure;imagesc(dff);axis image;
set(gca,'clim',[0,1]);

ss=size(dff_r);
merge=zeros([ss(1:2),3]);
merge(:,:,1)=dff_r;
merge(:,:,2)=dff_l;
merge(merge>1)=1;
merge(merge<0)=0;
figure;image(merge);


