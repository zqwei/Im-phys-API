% files=dir('Image*.tif');
% lastfile=files(end).name;
% ntrial=str2num(lastfile(end-6:end-4))-1;
% prefix=lastfile(1:end-7);
% 
% nimage=zeros(1,ntrial);
% for i=1:ntrial
%     name=[prefix,num2str(i,'%03d'),'.tif'];
%     a=dir(name);
%     if isempty(a)
%         nimage(i)=0;       
%     else
%         info=imfinfo(a.name);
%         nimage(i)=length(info);
%     end        
% end


%%
file=dir('*para.mat');
load(file(1).name,'ROI_list');
ROI_list=make_neuropil_mask(ROI_list);

for k=1:length(ROI_list)
    ROI_list(k).fmean=[];
end

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
% offset=121;
bev_range=(1:length(nimage))+offset;
trial=alltrial(bev_range);
save(file(1).name,'ROI_list','trial','nimage','offset','-append');

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


%%

frameRate=29.68/2;   %15HZ
range=-45:45;

% frameRate=29.68/4;   
% range=-20:30;



avg_l=zeros(512,512,length(range));
avg_r=zeros(512,512,length(range));

hdr=[];
f_fullfield=[];
rcount=0;
nimage1=[];
for i=1:ntrial
    name=[prefix,num2str(i,'%03d'),'.tif'];
    a=dir(name);
    if ~isempty(a)        
        [im h] = load_image(name);
        nimage1(i)=size(im,3);
        hdr=[hdr,h];
        cueframe=floor(frameRate*[trial(i).cuetime]);
        
        if (trial(i).type=='r')&&(trial(i).correct)&&~(trial(i).early)            
            trial_im=im(:,:,cueframe+range);            
            avg_r=avg_r+trial_im;
            rcount=rcount+1;
        elseif (trial(i).type=='l')&&(trial(i).correct)&&~(trial(i).early)            
            trial_im=im(:,:,cueframe+range);            
            avg_l=avg_l+trial_im;
        end
        
        im=reshape(im,[],nimage(i));
        tic
        for k=1:length(ROI_list)
            ROI_list(k).fmean=[ROI_list(k).fmean;mean(im(ROI_list(k).pixel_list,:))'];
            ROI_list(k).fmean_neuropil=[ROI_list(k).fmean_neuropil;mean(im(ROI_list(k).neuropil_list,:))'];
        end
        f_fullfield=[f_fullfield;mean(im)'];
    end
%     for k=1:length(ROI_list_PT)
%         ROI_list_PT(k).fmean=[ROI_list_PT(k).fmean;mean(im(ROI_list_PT(k).pixel_list,:))'];        
%     end
    toc
    i
end


%% merge movie

base=1:5;
bin_r=avg_r(1:2:end,1:2:end,:)+avg_r(1:2:end,2:2:end,:)+avg_r(2:2:end,1:2:end,:)+avg_r(2:2:end,2:2:end,:);
bin_l=avg_l(1:2:end,1:2:end,:)+avg_l(1:2:end,2:2:end,:)+avg_l(2:2:end,1:2:end,:)+avg_l(2:2:end,2:2:end,:);

f0_r=mean(bin_r(:,:,base),3)+1;
f0_l=mean(bin_l(:,:,base),3)+1;

f0_r=repmat(f0_r,[1,1,size(avg_r,3)]);
f0_l=repmat(f0_l,[1,1,size(avg_r,3)]);
dff_r=(bin_r-f0_r)./f0_r;
dff_l=(bin_l-f0_l)./f0_l;
ss=size(dff_r);
merge=zeros([ss,3]);
merge(:,:,:,1)=dff_r;
merge(:,:,:,2)=dff_l;
merge(merge>2)=2;
merge(merge<-2)=-2;
merge=merge/2;

%%
save(file(1).name,'ROI_list','f_fullfield','avg_r','avg_l','-append');