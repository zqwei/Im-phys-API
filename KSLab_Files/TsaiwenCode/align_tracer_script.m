prefix='an026_2013_08_12';
depth=465;


im=readTifStack('session_pertrial_mean_chan_01.tif');
im_trialmean=im;
avg=mean(im_trialmean,3);
save([prefix,'_',num2str(depth),'_para.mat'],'depth','avg','im_trialmean');
%%
file=dir('*_red*.tif');
if ~isempty(file)
    zstack=readTifStack(file(1).name);
    ss=size(zstack);
    zstack=reshape(zstack,ss(1),ss(2),2,[]);        
    green=squeeze(zstack(:,:,1,:));
    red=squeeze(zstack(:,:,2,:));
    red=red-min(red(:));
    srange=20;
    mcorr=[];
    for i=3:12%1:size(green,3)
        [Ireg,shift1,shift2,maxcorr]=imshiftreg(green(:,:,i),avg,srange);
        mcorr(i)=maxcorr;
    end
    [M,idx]=max(mcorr)
%     idx=7;
    [Ireg,shift1,shift2,maxcorr,xcorr]=imshiftreg(green(:,:,idx),avg,srange);
    
    g=mean(green(:,:,idx+(-1:1)),3);
    r=mean(red(:,:,idx+(-2:2)),3);
    green_reg=imshift(g,shift1,shift2);
    red_reg=imshift(r,shift1,shift2);
    id=green_reg>0;
    
    figure;plot(green_reg(id),red_reg(id),'.');
    b = robustfit(green_reg(id),red_reg(id));
    hold on;
    x=[min(green_reg(id)),max(green_reg(id))];
    plot(x,b(1)+x*b(2),'r');
    red_reg(id)=red_reg(id)-b(2)*green_reg(id);
    
    overlay=zeros(ss(1),ss(2),3);
    MM=prctile(green_reg(:),99.9);
    mm=mode(green_reg(:));
    overlay(:,:,2)=(green_reg-mm)/(MM-mm);
    MM=prctile(red_reg(:),99.5);
%     mm=mode(red_reg(:))*0.8;
    mm=prctile(red_reg(:),12);
    overlay(:,:,1)=(red_reg-mm)/(MM-mm);
    overlay(overlay>1)=1;
    overlay(overlay<0)=0;
    imcompare2(avg,overlay);
    overlay_red=overlay;
    zstack_red=zstack;
    save([prefix,'_',num2str(depth),'_para.mat'],'zstack_red','overlay_red','-append')
end
%%
file=dir('*_farred*.tif');
if ~isempty(file)
    zstack=readTifStack(file(1).name);
    ss=size(zstack);
    zstack=reshape(zstack,ss(1),ss(2),2,[]);        
    green=squeeze(zstack(:,:,1,:));
    red=squeeze(zstack(:,:,2,:));
    red=red-min(red(:));
    srange=20;
    mcorr=[];
    for i=3:12%1:size(green,3)
        [Ireg,shift1,shift2,maxcorr]=imshiftreg(green(:,:,i),avg,srange);
        mcorr(i)=maxcorr;
    end
    [M,idx]=max(mcorr)
    [Ireg,shift1,shift2,maxcorr,xcorr]=imshiftreg(green(:,:,idx),avg,srange);
    
    g=mean(green(:,:,idx+(-1:1)),3);
    r=mean(red(:,:,idx+(-2:2)),3);
    green_reg=imshift(g,shift1,shift2);
    red_reg=imshift(r,shift1,shift2);
    
    overlay=zeros(ss(1),ss(2),3);
    MM=prctile(green_reg(:),99.9);
    mm=mode(green_reg(:));
    overlay(:,:,2)=(green_reg-mm)/(MM-mm);
    MM=prctile(red_reg(:),99.5);
%     mm=mode(red_reg(:))*0.8;
    mm=prctile(red_reg(:),12);
    overlay(:,:,1)=(red_reg-mm)/(MM-mm);
    overlay(overlay>1)=1;
    overlay(overlay<0)=0;
    imcompare2(avg,overlay);
    overlay_farred=overlay;
    zstack_farred=zstack;
    save([prefix,'_',num2str(depth),'_para.mat'],'overlay_farred','zstack_farred','-append');
end


% 
% %%
% load('ROI_list');
% rval=[];
% gval=[];
% for i=1:length(ROI_list)
%     pixel_list=ROI_list(i).pixel_list;
%     rval(i)=mean(red_reg(pixel_list));
%     gval(i)=mean(green_reg(pixel_list));
% end
% 
% %%
% b = robustfit(gval,rval);
% figure;plot(gval,rval,'.');
% hold on;
% x=[0,max(gval)];
% plot(x,b(1)+b(2)*x,'-');
% rval_pre=b(1)+b(2)*gval;
% rval_comp=rval-rval_pre;