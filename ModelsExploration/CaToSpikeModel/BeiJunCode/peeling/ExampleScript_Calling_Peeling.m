function [out] = Deconv_Peeling(out)
%
% Peeling algrithm to extract spike trains and instantanuoes spike rates underlying calcium
% transients; 
% Output: spike train vectors and inst. rate vectors for all cells and all
% trials
% last update 8.2.2104, Fritjof Helmchen, Brain Research Institute, UZH
%
frameRate = out{1,1}.Ca.sample_rate;
convert2Hz = frameRate;  % conversion factor to instantaneous rate

% presumed elementary calcium transient, onset exponential and single-exponential decay
% values for YC-Nano140
drr_amp = 0.045;  % DR/R amplitude 4.5% 
drr_onset = 0.05;  % in sec
drr_tau = 0.7;  % in sec

cellNum = numel(out{1,1}.Ca.roiLabel);
trialNum = numel(out);
totNum = cellNum*trialNum;

allspike_out = cell(cellNum,1);
allspike_density = cell(cellNum,1);
    
h = waitbar(0,'Peeling! Please wait...');

for ii = 1:trialNum
    orgFrames = numel(out{1,ii}.Motor.Params_18(:,2));
    extraFrames = frameRate*2; % equivalent to 2 sec
    totFrames = orgFrames  + extraFrames;

    % presumed elementary calcium transient
    kern = drr_amp*(1-exp(-((1:totFrames)*1-1)/frameRate/drr_onset)).*exp(-((1:totFrames)*1-1)/frameRate/drr_tau);

    % prepare arrays 
    Ca_tmp = zeros(1,totFrames); % expanded Ca transient
    spike_out = zeros(1,orgFrames);
    spike_density = zeros(1,orgFrames);
    
    % smoothing kernel for converting spike trains into inst. firing rate (in Hz);
    % sigma in units of frames (default: 1 frames) 
    gausskern = normpdf((1:orgFrames),round(0.5*orgFrames),1);
    
    for numi= 1:cellNum
        Ca_tmp(1,1:orgFrames) = out{1,ii}.Ca.dRR{numi,1}(1,1:orgFrames);
        mean_last = sum(out{1,ii}.Ca.dRR{numi,1}(1,(orgFrames-2):orgFrames))/3;
        Ca_tmp(1,(orgFrames+1):(orgFrames+extraFrames)) = mean_last*exp(-((1:extraFrames)*1-1)/frameRate/drr_tau);
    
        inan = find(isnan(Ca_tmp(1,:)));
        Ca_tmp(1,inan) = 0;  % remove potentials NaN entrys and replace with zero
        
        % apply Peeling algorithm to expanded Ca transient
        pdrr = Ca_tmp(1,1:orgFrames)*100;
        
        [ca_p, exp_p, peel_p, data] = InitPeeling(pdrr, frameRate);
        [ca_p, peel_p, data] = Peeling(pdrr, frameRate);
        
        spike_out(1,1:orgFrames) = data.spiketrain';
        APdens = conv(data.spiketrain', gausskern);
        spike_density(1,1:orgFrames) = convert2Hz*APdens(round(0.5*orgFrames):round(0.5*orgFrames)+orgFrames-1);
        
        allspike_out{numi,1}= spike_out;
        allspike_density{numi,1} = spike_density;
        
        waitbar(((ii-1)*cellNum+numi)/totNum);
    end
    out{1,ii}.Ca.Peel_Spikes = allspike_out;
    out{1,ii}.Ca.Peel_InstRate = allspike_density;
end

close(h); 








