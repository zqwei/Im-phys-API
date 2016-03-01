function S = modelCalcium(S,doPlot)
% see ParseInputs function for input arguments and default parameters
% S ... configuration structure
% doPlot ... generate plots?
% run the function without any input arguments first, to get the configuration structure S:
% S = modelCalcium
% parameters in S can then be changed (see ParseInputs function for explanation of the parameters)

% Henry Luetcke & Fritjof Helmchen
% Brain Research Institut
% University of Zurich
% Switzerland
%
% For details on the simulations please also see the corresponding manuscript:
% Luetcke, Gerhard, Zenke, Gerstner & Helmchen (in revision), Frontiers in Neural Circuits

%% Parameters
if ~nargin
    S = struct;
    doPlot = 1;
end

if ~isstruct(S)
    error('Input must be either structure or string variable')
end

[S,cellNo] = ParseInputs(S);

cmap = lines(cellNo);
if cellNo == 1
    cmap = [0 0 0];
end

% simulation should be 1 s or so longer (allow for steady-state in case of high
% firing rates)
dur = S.dur + S.offset;

x = 1/S.samplingRate:1/S.samplingRate:dur;
titleStr = 'Simulation';
if doPlot
    fig = figure('Name',titleStr,'NumberTitle','off','Color','white');
    hold all
end
ca = zeros(cellNo,numel(x));
dff = zeros(cellNo,numel(x));
spkTCell = S.spkTimes;

ca_genmode = S.ca_genmode;
spk_recmode = S.spk_recmode; 
tauOn = S.tauOn;
A1 = S.A1;
tau1 = S.tau1;
A2 = S.A2;
tau2 = S.tau2;

% expected peak amplitude based on exponential rise and decay
switch lower(ca_genmode)
    case 'lindff'
        PeakA = A1 .* (tau1./tauOn.*(tau1./tauOn+1).^-(tauOn./tau1+1));
    case 'satdff'
        S.ca_amp1=S.ca_amp./(1+S.ca_kappas+S.kappab);             % set for consistency
        S.ca_tau1=(1+S.ca_kappas+S.kappab)./S.ca_gamma;
        PeakCa = S.ca_amp1 + S.ca_rest;
        PeakDFF = Calcium2Fluor(PeakCa,S.ca_rest,S.kd,S.dffmax);
        PeakA = PeakDFF .* ...
            (S.ca_tau1./S.ca_onsettau.*(S.ca_tau1./S.ca_onsettau+1).^-(S.ca_onsettau./S.ca_tau1+1));
    otherwise
        error('Calcium trace generation mode illdefined');
end
sdnoise = PeakA./S.snr; % determines how much Gaussian noise is added to the DFF trace

samplingRate = S.samplingRate;
offset = S.offset;
spkTCell2 = cell(1,cellNo);

% in principle it is possible to run simulations for multiple neurons at the same time
% in this case input arguments that differ from cell to cell have to be provided in vectorized
% format (e.g. one entry per cell)
% this is only supported for parameters where this makes experimentally sense, e.g. for the firing
% rate, decay time etc. but not for imaging parameters, such as the frame rate
for m = 1:cellNo
    currentSpkT = spkTCell{m};
    if isempty(spkTCell{m})
        fprintf('No spikes for cell %1.0f. Skipping ...\n',m);
        if m == 1
            close all
            error('Must have spikes for cell 1. Try to increase firing rate / simulation duration.');
        end
        continue
    end
    if isfield(S,'A1sigma') && ~isempty(S.A1sigma) && isfield(S,'tau1sigma') && ~isempty(S.tau1sigma)
        % convolution for each spike (slower, allows variable model calcium transient)
        for n = 1:numel(currentSpkT)
            currentA1 = random('normal',A1(m),A1(m).*S.A1sigma);
            currentTau1 = random('normal',tau1(m),tau1(m).*S.tau1sigma);
            y = spkTimes2Calcium(currentSpkT(n),tauOn(m),currentA1,currentTau1,A2(m),...
                tau2(m),samplingRate,dur);
            dff(m,:) = dff(m,:) + y(1:numel(x));
        end
    else
        switch lower(ca_genmode)
            case 'lindff'
                % convolution over all spikes (faster, same model calcium transient)
                modelTransient = spkTimes2Calcium(0,tauOn(m),A1(m),tau1(m),A2(m),...
                    tau2(m),samplingRate,dur);
                spkVector = zeros(1,numel(x));
                for i = 1:numel(currentSpkT)
                    [~,idx] = min(abs(currentSpkT(i)-x));
                    spkVector(idx) = spkVector(idx)+1;
                end
                dffConv = conv(spkVector,modelTransient);
                dff(m,:) = dffConv(1:length(x));
            case 'satdff'
                % taking saturation into account by solving the single comp model differential equation
                % piecewise, then applying nonlinear transformation from ca to dff
                ca = spkTimes2FreeCalcium(currentSpkT,S.ca_amp,S.ca_gamma,S.ca_onsettau,S.ca_rest, S.ca_kappas,...
                    S.kd, S.conc,S.samplingRate,dur);
                dff(m,:) = Calcium2Fluor(ca,S.ca_rest,S.kd,S.dffmax);
        end
    end
    currentSpkT = currentSpkT - offset;
    currentSpkT(currentSpkT<0) = [];
    spkTCell2{m} = currentSpkT;
end
S.spkTimes = spkTCell2;

subPlotNo = length(S.frameRate)+2;
axisHandles = [];

%% DFF
xPlot = x - S.offset;
dff(:,xPlot<0) = [];
xPlot(xPlot<0) = [];
if doPlot
    subplot(subPlotNo,1,1)
    yOffset =0;
    for m = 1:cellNo
        if isempty(spkTCell{m})
            continue
        end
        if m > 1
            currentDff = dff(m,:)+(yOffset-min(dff(m,:)));
        else
            currentDff = dff(m,:);
        end
        plot(xPlot,currentDff,'Color',cmap(m,:)); hold on
        yOffset = max(currentDff)+max(dff(:))/5;
        for n = 1:numel(S.spkTimes{m})
            h_err = errorbar(S.spkTimes{m}(n),max(currentDff),0,...
                max(dff(:))/10,'Color',cmap(m,:));
            removeErrorBarEnds(h_err)
            set(h_err,'LineWidth',2)
        end
    end
    set(gca,'ylim',[min(dff(:)) max(currentDff)+max(dff(:))/10],'box','off')
    axisHandles(1) = gca;
    ylabel('DFF / %');
end
S.data.dff = dff;
S.data.ca = ca;

%% noisy DFF
noisyDFF = zeros(cellNo,size(dff,2));
for m = 1:cellNo
    if isempty(spkTCell{m})
        continue
    end
    
    whiteNoise = sdnoise(m).*randn(1,size(dff,2));
    % it is possible to 'color' the noise, to mimick possible temporal autocorrelations (not
    % extensively tested)
    %     a = 0.99; % 0.99 gives auto-cov similar to exp. data
    % First Order Low pass filter y(n)=a*y(n-1)+(1-a)*x(n)
    % Filter Transfer function Y(Z) = X(Z)*(1-a)/(1-aZ^-1)
    %     coloredNoise = (filter(1-a,[1 -a],whiteNoise)) + whiteNoise;
    coloredNoise = whiteNoise;
    noisyDFF(m,:) = dff(m,:) + coloredNoise;
end
if doPlot
    subplot(subPlotNo,1,2)
    yOffset =0;
    for m = 1:cellNo
        if isempty(spkTCell{m})
            continue
        end
        if m > 1
            currentDff = noisyDFF(m,:)+(yOffset-min(noisyDFF(m,:)));
        else
            currentDff = noisyDFF(m,:);
        end
        plot(xPlot,currentDff,'Color',cmap(m,:)); hold on
        yOffset = max(currentDff)+max(noisyDFF(:))/5;
        for n = 1:numel(S.spkTimes{m})
            h_err = errorbar(S.spkTimes{m}(n),max(currentDff),0,...
                max(noisyDFF(:))/10,'Color',cmap(m,:));
            removeErrorBarEnds(h_err)
            set(h_err,'LineWidth',2)
        end
    end
    set(gca,'ylim',[min(noisyDFF(:)) ...
        max(currentDff)+max(noisyDFF(:))/10],'box','off')
    axisHandles(2) = gca;
    ylabel('DFF / %');
end
S.data.noisyDFF = noisyDFF;

%% noisy DFF at frame rate
lowResT = 1/S.frameRate:1/S.frameRate:max(xPlot);
lowResT = (1/S.frameRate:1/S.frameRate:max(xPlot)) - (0.5 .* 1/S.frameRate);
[~,idxList] = findClosest(lowResT,xPlot);
noisyDFFlowRes = zeros(m,numel(idxList));
for m = 1:cellNo
    if isempty(spkTCell{m})
        continue
    end
    noisyDFFlowRes(m,:) = noisyDFF(m,idxList);
end

if doPlot
    subplot(subPlotNo,1,3)
    yOffset = 0;
    for m = 1:cellNo
        if isempty(spkTCell{m})
            continue
        end
        if m > 1
            currentDff = noisyDFFlowRes(m,:)+(yOffset-min(noisyDFF(m,:)));
        else
            currentDff = noisyDFFlowRes(m,:);
        end
        plot(lowResT,currentDff,'Color',cmap(m,:)); hold on
        yOffset = max(currentDff)+max(noisyDFFlowRes(:))/5;
        for n = 1:numel(S.spkTimes{m})
            h_err = errorbar(S.spkTimes{m}(n),max(currentDff),0,...
                max(noisyDFFlowRes(:))/10,'Color',cmap(m,:));
            removeErrorBarEnds(h_err)
            set(h_err,'LineWidth',2)
        end
    end
    set(gca,'ylim',[min(noisyDFFlowRes(:)) ...
        max(currentDff)+max(noisyDFFlowRes(:))/10],'box','off')
    axisHandles(3) = gca;
    ylabel(sprintf('DFF (%1.0f Hz)',S.frameRate));
    xlabel('Time / s');
    linkaxes(axisHandles,'x')
end
dff = zeros(cellNo,length(noisyDFFlowRes));
for m = 1:cellNo
    if isempty(spkTCell{m})
        continue
    end
    % now perform baseline setting based on ratio of expected and real SD
    % of trace
    ratioSD = std(noisyDFFlowRes(m,:))./sdnoise(m);
    
    % empirical function describing baseperc as function of ratioSD
    basePerc = 50 - 50 / (1+exp(-(ratioSD-1-S.snr(m)*0.25)*10/S.snr(m)));

    baseValue = prctile(noisyDFFlowRes(m,:),basePerc);
    % sometimes it is better to specify the baseline explicitely, for example in cases of short
    % bursts of activity
    baseValue = 0;
    dff(m,:) = noisyDFFlowRes(m,:) - baseValue;
end

S.data.noisyDFFlowResT = noisyDFFlowRes;

if doPlot
    titleStr = sprintf('Reconstruction - %s',S.reconAlg);
    fig = figure('Name',titleStr,'NumberTitle','off','Color','white'); hold all
end


%% Reconstruction
peelOpts = S.peelOptions;
% schmitt = [1.75 -1 0.3];
S.recon = struct;
switch lower(S.reconAlg)
    case 'peeling'
        if isfield(S,'schmitt')
            schmitt = S.schmitt;
        end
        [S.recon.spikePredict,S.recon.peel] = ...
            doPeel(dff,S.frameRate,S.ca_genmode, S.spk_recmode,S.A1,S.tau1,S.A2,S.tau2,S.tauOn,S.ca_rest, S.ca_amp, S.ca_gamma, S.ca_onsettau, ...
                   S.ca_amp1,S.ca_tau1,S.ca_kappas, S.dffmax, S.kd, S.kappab, sdnoise, peelOpts.schmitt, peelOpts.peel_p);
end

dffRecon = zeros(cellNo, length(lowResT));
S.recon.spkTimesRecon = cell(1,cellNo);
for m = 1:cellNo
    if isempty(spkTCell{m})
        continue
    end
    S.recon.spkTimesRecon{m} = S.recon.spikePredict{m};
    
    % reconstructed model trace based on reconstructed spike times
    switch lower(S.spk_recmode)
        case 'lindff'
            modelTransient = spkTimes2Calcium(0,tauOn(m),A1(m),tau1(m),A2(m),...
                tau2(m),S.frameRate,dur);
            spkVector = zeros(1,numel(lowResT));
            for i = 1:numel(S.recon.spkTimesRecon{m})
                [~,idx] = min(abs(S.recon.spkTimesRecon{m}(i)-lowResT));
                spkVector(idx) = spkVector(idx)+1;
            end
            dffConv = conv(spkVector,modelTransient);
            dffRecon(m,:) = dffConv(1:length(lowResT));
        case 'satdff'
            ca = spkTimes2FreeCalcium(S.recon.spkTimesRecon{m},S.ca_amp,S.ca_gamma,S.ca_onsettau,S.ca_rest, S.ca_kappas,...
                S.kd, S.conc,S.frameRate,dur);
            dfftmp(m,:) = Calcium2Fluor(ca,S.ca_rest,S.kd,S.dffmax);
            dffRecon(m,:) = dfftmp(1:length(lowResT));
        otherwise
            error('Model trace generation failed. Undefined SpkReconMode.');
    end
end

if doPlot
    yOffset =0;
    for m = 1:cellNo
        if isempty(spkTCell{m})
            continue
        end
        if m > 1
            currentDff = dff(m,:)+(yOffset-min(dff(m,:)));
        else
            currentDff = dff(m,:);
        end
        plot(lowResT,currentDff,'Color',[0 0 0]); hold on       
        plot(lowResT,dffRecon,'Color',[0 1 0]);
        yOffset = max(currentDff)+max(dff(:))/5;
        for n = 1:numel(S.spkTimes{m})
            h_err = errorbar(S.spkTimes{m}(n),max(currentDff),0,...
                max(dff(:))/10,'Color',[0 0 0]);
            removeErrorBarEnds(h_err)
            set(h_err,'LineWidth',2)
        end
        for n = 1:numel(S.recon.spkTimesRecon{m})
            h_err = errorbar(S.recon.spkTimesRecon{m}(n),...
                max(currentDff)+max(dff(:))/10,0,max(dff(:))/10,...
                'Color',[1 0 0]);
            removeErrorBarEnds(h_err)
            set(h_err,'LineWidth',2)
        end
    end
    set(gca,'ylim',...
        [min(dff(:)) max(currentDff)+2*max(dff(:))/10],'box','off')
    ylabel('DFF / %s');
    xlabel('Time / s');
end

end

%% Function - doPeel
function [spikePredict,peel] = doPeel(dff,frameRate,ca_genmode,spk_recmode, A1,tau1,A2,tau2,tauOn,ca_rest, ca_amp,ca_gamma, ca_onsettau, ca_amp1,ca_tau1,...
    ca_kappas, dffmax, kd, kappab, sdnoise, schmitt, peel_p)
for m = 1:size(dff,1)
    %     [ca_p,exp_p,peel_p, data] = InitPeeling(dff(m,:),frameRate);
    
    % override some parameters
    ca_p.ca_genmode = ca_genmode;
    ca_p.amp1=A1(m);
    ca_p.tau1=tau1(m);
    ca_p.amp2=A2(m);
    ca_p.tau2=tau2(m);
    ca_p.onsettau = tauOn(m);
    ca_p.ca_rest= ca_rest(m);
    ca_p.ca_amp= ca_amp(m);
    ca_p.ca_gamma= ca_gamma(m);
    ca_p.ca_onsettau= ca_onsettau(m);
    ca_p.ca_amp1= ca_amp1(m);
    ca_p.ca_tau1= ca_tau1(m);
    ca_p.ca_kappas= ca_kappas(m);

    exp_p.dffmax=dffmax(m);
    exp_p.kd=kd(m);
    exp_p.kappab=kappab(m);
    exp_p.noiseSD = sdnoise(m);
    
    switch lower(ca_p.ca_genmode)
        case 'lindff'
            PeakA = A1 .* (tau1./tauOn.*(tau1./tauOn+1).^-(tauOn./tau1+1));
        case 'satdff'
            ca_p.ca_amp1=ca_p.ca_amp./(1+ca_p.ca_kappas+exp_p.kappab);             % set for consistency
            ca_p.ca_tau1=(1+ca_p.ca_kappas+exp_p.kappab)./ca_p.ca_gamma;
            PeakCa = ca_p.ca_amp1 + ca_p.ca_rest;
            PeakDFF = Calcium2Fluor(PeakCa,ca_p.ca_rest,exp_p.kd,exp_p.dffmax);
            PeakA = PeakDFF .* ...
                (ca_p.ca_tau1./ca_p.ca_onsettau.*(ca_p.ca_tau1./ca_p.ca_onsettau+1).^-(ca_p.ca_onsettau./ca_p.ca_tau1+1));
        otherwise
            error('Error in DoPeel. Illdefined SpkGenMode');
    end
    
    snr = PeakA./exp_p.noiseSD;
    
    peel_p.spk_recmode = spk_recmode;
    peel_p.smtthigh = schmitt(1)*exp_p.noiseSD;
    peel_p.smttlow = schmitt(2)*exp_p.noiseSD;
    
    peel_p.smttmindur= schmitt(3);
    peel_p.slidwinsiz = 10.0;
    peel_p.fitupdatetime=0.5;
    
    peel_p.doPlot = 0;
    
    [ca_p, peel_p, data] = Peeling(dff(m,:), frameRate, ca_p, exp_p, peel_p);
    
    spikePredict{m} = data.spikes;
    peel{m} = data.peel;
end

end

%% Function - ParseInputs
function [S, doVectorized] = ParseInputs(S)
doVectorized = 0;
% note that all vectorizable inputs must be the same length (or length=1)

% if ~isfield(S,'usefreeca')
%    S.usefreeca = 0; % flag, 0 - off, 1 - Ca calculation including dye saturation switched on
% end

if ~isfield(S,'ca_genmode')
    S.ca_genmode = 'linDFF'; %  default Ca trace generation mode is most simple linear case 
end

if ~isfield(S,'spk_recmode')
    S.spk_recmode = 'linDFF'; %  default spk reconstruction mode is most simple linear case 
end

if ~isfield(S,'ca_onsettau')
    S.ca_onsettau = 0.02; % Ca transient onset tau (s) / vectorizable for network of neurons
end
if length(S.ca_onsettau) > 1
    doVectorized = 1;
end

if ~isfield(S,'ca_amp')
    S.ca_amp = 7600; % single AP Ca transient amplitude (nM) / vectorizable for network of neurons
end
if length(S.ca_amp) > 1
    doVectorized = 1;
end

if ~isfield(S,'ca_gamma')
    S.ca_gamma = 400; % Ca extrusion rate (1/s) / vectorizable for network of neurons
end
if length(S.ca_gamma) > 1
    doVectorized = 1;
end

if ~isfield(S,'ca_amp1')
    S.ca_amp1 = 0; % Ca transient amp (nM), to be calculated / vectorizable for network of neurons
end
if length(S.ca_amp1) > 1
    doVectorized = 1;
end

if ~isfield(S,'ca_tau1')
    S.ca_tau1 = 0; % Ca transient decay time constant (s), to be calculated / vectorizable for network of neurons
end
if length(S.ca_tau1) > 1
    doVectorized = 1;
end

if ~isfield(S,'ca_kappas')
    S.ca_kappas = 100; % endogenous Ca-binding ratio
end

if ~isfield(S,'ca_rest')
    S.ca_rest = 50; % Ca restng level (nM) / 
end

if ~isfield(S,'dffmax')
    S.dffmax = 93; % saturating dff max (in percent) / 
end

if ~isfield(S,'kd')
    S.kd = 250; % dye dissociation constant (nM)
end

if ~isfield(S,'conc')
    S.conc = 50000; % dye total concentration (nM)
end

if ~isfield(S,'kappab')
    S.kappab = S.kd.*S.conc./(S.ca_rest+S.kd).^2;  % exogenous Ca-binding ratio
end

if (S.ca_genmode == 'linDFF')
elseif (S.ca_genmode == 'satDFF')
    S.ca_amp1=S.ca_amp./(1+S.ca_kappas+S.kappab);             % init for consistency
    S.ca_tau1=(1+S.ca_kappas+S.kappab)./S.ca_gamma;
end

if ~isfield(S,'A1')
    S.A1 = 8.5; % single AP DFF in % / vectorizable for network of neurons
end
if length(S.A1) > 1
    doVectorized = 1;
end

if ~isfield(S,'tau1')
    S.tau1 = 0.5; % indicator decay time in s / vectorizable for network of neurons
end
if length(S.tau1) > 1
    doVectorized = 1;
end

if ~isfield(S,'A1sigma')
    S.A1sigma = []; % spread of A1 around mean
end

if ~isfield(S,'tau1sigma')
    S.tau1sigma = []; % spread of tau1 around mean
end

if ~isfield(S,'A2')
    S.A2 = 0; % second amplitude for double-exp decay / vectorizable for network of neurons
end
if length(S.A2) > 1
    doVectorized = 1;
end

if ~isfield(S,'tau2')
    S.tau2 = 1.0; % decay time for second exponential / vectorizable for network of neurons
end
if length(S.tau2) > 1
    doVectorized = 1;
end

if ~isfield(S,'tauOn')
    S.tauOn = 0.01; % onset time in s / vectorizable for network of neurons
end
if length(S.tauOn) > 1
    doVectorized = 1;
end

if ~isfield(S,'dur')
    S.dur = 30; % simulation duration in s
end

if ~isfield(S,'spikeRate')
    S.spikeRate = 0.2; % firing rate in Hz / vectorizable for network of neurons
end
if length(S.spikeRate) > 1
    doVectorized = 1;
end

if ~isfield(S,'snr')
    S.snr = 5; % SNR (=A/SDnoise) / vectorizable for network of neurons
end
if length(S.snr) > 1
    doVectorized = 1;
end

if ~isfield(S,'samplingRate')
    S.samplingRate = 1000; % rate for calculation of calcium concentration and DFF in Hz
end
if ~isfield(S,'frameRate')
    S.frameRate = 30; % simulate imaging of DFF with this rate in Hz
end

if ~isfield(S,'offset')
    S.offset = 1; % in s
end
if ~isfield(S,'maxdt')
    S.maxdt = 0.5; % for evaluation, in s
end
if ~isfield(S,'reconAlg')
    S.reconAlg = 'peeling'; % reconstruction algorithm ('peeling' or 'none')
end

if ~isfield(S,'peelOptions')
    peelOptions.schmitt = [1.75 -1 0.3];
    peelOptions.peel_p.optimizeSpikeTimes = 1;
    peelOptions.peel_p.fitonset = 0;
    S.peelOptions = peelOptions;
end

switch S.reconAlg
    case 'none'
        warning('Skipping reconstruction (and evaluation).');
end

if ~isfield(S,'recycleSpikeTimes')
    S.recycleSpikeTimes = 0; % reuse existing spike train
end

% if ~isfield(S,'doEvaluation')
%     S.doEvaluation = 1; % compare reconstructed and original spike trains
% end


% check if vectorized arguments are consistent
if doVectorized
    argVector = {S.A1 S.tau1 S.A2 S.tau2 S.tauOn S.spikeRate S.snr};
    if length(unique(cellfun(@length,argVector))) > 2 % length can be either 1 or number of neurons
        error('Argument lengths are inconsistent');
    end
    doVectorized = unique(cellfun(@length,argVector));
    if doVectorized(1) == 1
        doVectorized = doVectorized(2);
    else
        doVectorized = doVectorized(1);
    end
    if length(S.A1) == 1
        S.A1 = repmat(S.A1,1,doVectorized);
    end
    if length(S.tau1) == 1
        S.tau1 = repmat(S.tau1,1,doVectorized);
    end
    if length(S.A2) == 1
        S.A2 = repmat(S.A2,1,doVectorized);
    end
    if length(S.tau2) == 1
        S.tau2 = repmat(S.tau2,1,doVectorized);
    end
    if length(S.tauOn) == 1
        S.tauOn = repmat(S.tauOn,1,doVectorized);
    end
    if length(S.spikeRate) == 1
        S.spikeRate = repmat(S.spikeRate,1,doVectorized);
    end
    if length(S.snr) == 1
        S.snr = repmat(S.snr,1,doVectorized);
    end
else
    doVectorized = 1;
end

if ~isfield(S,'spkTimes') || ~S.recycleSpikeTimes
    S.spkTimes = cell(doVectorized,1);
    for n = 1:doVectorized
        %     poisson spike train
        S.spkTimes{n,1} = PoissonSpikeTrain(S.spikeRate(n),S.dur+S.offset);
    end
else
    if ~iscell(S.spkTimes)
        spkTimes = S.spkTimes;
        S.spkTimes = cell(1,1);
        S.spkTimes{1} = spkTimes;
    end
    if length(S.spkTimes) ~= doVectorized
        warning('Spike time vector length not consistent with vectorized arguments. Setting to spike time vector length.');
        doVectorized = length(S.spkTimes);
        S.A1 = repmat(S.A1(1),1,doVectorized);
        S.tau1 = repmat(S.tau1(1),1,doVectorized);
        S.A2 = repmat(S.A2(1),1,doVectorized);
        S.tau2 = repmat(S.tau2(1),1,doVectorized);
        S.tauOn = repmat(S.tauOn(1),1,doVectorized);
        S.snr = repmat(S.snr(1),1,doVectorized);
        % set S.spikeRate to average firing rate
        for n = 1:doVectorized
            avgFiringRate = numel(S.spkTimes{n}) / S.dur;
            S.spikeRate(n) = avgFiringRate;
        end
    end
    for n = 1:doVectorized
        S.spkTimes{n} = S.spkTimes{n} + S.offset;
    end
end

end


