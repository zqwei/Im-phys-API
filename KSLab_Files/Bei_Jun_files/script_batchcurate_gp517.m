files=dir('*para.mat');

for i=1:length(files)
    load(files(i).name,'para');
    if ~isfield(para,'t_frame')
        para.t_frame=(1:(para.recording_len*para.fs1))'/para.fs1;
        para.t_ephys=(1:(para.recording_len*para.fs2))'/para.fs2;
    end
    
%     dff=get_baseline_corr_dff(para);
    f0=mode(para.fmean);
    % f0=30;
    dff=(para.fmean-f0)/f0;
    spk{1,1}=para.t_ephys(para.peakcurate>0);

   
    %% one step method
    para_start=[20   20    1.8    0.5   0.02];
    vector01=fcn_perispike(para.peakcurate,dff,para.t_ephys,para.t_frame);
    [para_final,mse_final]=gcamp6_model_4para(spk,dff,vector01,para.t_frame,para_start);
   
    
    firing_rate=sum(para.peakcurate)/para.recording_len;
    
    
    [CaTraces]=func_getCaTraces_general_new(spk,para.t_frame,para_final);
    
    %vector01=fcn_perispike(para.peak,dff,para.t_ephys,para.t_frame);
    dff_mod=dff(vector01);
    CaTraces_mod=CaTraces(vector01);
    coreff=corrcoef(dff_mod,CaTraces_mod);
    coreff=coreff(1,2);
    fiterr=sum(((dff_mod-CaTraces_mod).^2))./sum((dff_mod.^2));

    
%     corr_model_data=corrcoef(dff,CaTraces);
%     corr_model_data=corr_model_data(1,2);
%     fittingidx=sum(((dff-CaTraces).*(dff-CaTraces)))./sum((dff.*dff));
    fitresult=[dff,CaTraces];
    figure('unit','inches','position',[3,3,8,2.8]);hold on;
    plot(para.t_frame,dff,'b');plot(para.t_frame,CaTraces,'r');
    plot(para.t_frame,dff-CaTraces-3,'k')
    plot(para.t_ephys,para.filt*0.8-6,'color',[1,1,1]*0.5,'linewidth',0.5);
    spk_time=find(para.peakcurate);
    plot(para.t_ephys(spk_time),-5,'rx');
    title(['corrcoef = ',num2str(coreff,'%0.3f'),' t_rise = ',num2str(para_final(5),'%0.3f'),' t_decay = ',num2str(para_final(4),'%0.3f')]);
     set(gca,'linewidth',2);
    set(gca,'fontsize',16);
    xlim([0,para.t_frame(end)]);
    %save(files(i).name,'para_final','firing_rate','coreff','fiterr','fitresult','-append');
end

