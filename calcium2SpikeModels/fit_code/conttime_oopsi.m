function SAMPLES = conttime_oopsi(Y,params)
    % Continuous time sampler
    % Y                     data (normalized in [0,1])
    % P                     intialization parameters (discrete time constant P.g required)
    % params                parameters structure
    % params.g              discrete time constant(s) (estimated if not provided)
    % params.sn             initializer for noise (estimated if not provided)
    % params.b              initializer for baseline (estimated if not provided)
    % params.c1             initializer for initial concentration (estimated if not provided)
    % params.Nsamples       number of samples after burn in (default 500)
    % params.B              number of burn in samples (default 200)
    % params.marg           flag for marginalized sampler (default 0)
    % params.upd_gam        flag for updating gamma (default 0)
    % params.gam_step       number of samples after which gamma is updated (default 50)
    % params.std_move       standard deviation of shifting kernel (default 3*Dt)
    % params.add_move       number of add moves per iteration (default T/100)
    % params.init           initial sample 
    % params.f              imaging rate (default 1)
    % params.p              order of AR model (p == 1 or p == 2, default 1)
    % params.defg           default discrete time constants in case constrained_foopsi cannot find stable estimates
    % params.TauStd         standard deviation for time constants in continuous time (default [0.2,2])
    % output struct SAMPLES
    % spikes                T x Nsamples matrix with spikes samples
    % bp                    Nsamples x 1 vector with samples for spiking prior probability
    % Am                    Nsamples x 1 vector with samples for spike amplitude
    % If marginalized sampler is used
    % Cb                    posterior mean and sd for baseline
    % Cin                   posterior mean and sd for initial condition
    % else
    % Cb                    Nsamples x 1 vector with samples for baseline
    % Cin                   Nsamples x 1 vector with samples for initial concentration
    % sn                    Nsamples x 1 vector with samples for noise variance
    % If gamma is updated
    % g                     Nsamples x p vector with the gamma updates
    % Author: Eftychios A. Pnevmatikakis and Josh Merel
    Y = Y(:);
    T = length(Y);
    % define default parameters
    defparams.g = [];
    defparams.sn = [];
    defparams.b = [];
    defparams.c1 = [];
    defparams.Nsamples = 400;
    defparams.B = 200;
    defparams.marg = 0;
    defparams.upd_gam = 1; 
    defparams.gam_step = 1;
    defparams.A_lb = 0.1*range(Y);
    defparams.b_lb = quantile(Y,0.01);
    defparams.c1_lb = 0;
    defparams.std_move = 3;
    defparams.add_move = ceil(T/100);
    defparams.init = [];
    defparams.f = 1;
    defparams.p = 1;
    defparams.defg = [0.6,0.95];
    defparams.TauStd = [.1,1];
    defparams.print_flag = 0;
    if nargin < 2
        params = defparams;
    else
        if ~isfield(params,'g'); params.g = defparams.g; end
        if ~isfield(params,'sn'); params.sn = defparams.sn; end
        if ~isfield(params,'b'); params.b = defparams.b; end
        if ~isfield(params,'c1'); params.c1 = defparams.c1; end
        if ~isfield(params,'Nsamples'); params.Nsamples = defparams.Nsamples; end
        if ~isfield(params,'B'); params.B = defparams.B; end
        if ~isfield(params,'marg'); params.marg = defparams.marg; end
        if ~isfield(params,'upd_gam'); params.upd_gam = defparams.upd_gam; end
        if ~isfield(params,'gam_step'); params.gam_step = defparams.gam_step; end
        if ~isfield(params,'std_move'); params.std_move = defparams.std_move; end
        if ~isfield(params,'add_move'); params.add_move = defparams.add_move; end
        if ~isfield(params,'init'); params.init = defparams.init; end
        if ~isfield(params,'f'); params.f = defparams.f; end
        if ~isfield(params,'p'); params.p = defparams.p; end
        if ~isfield(params,'defg'); params.defg = defparams.defg; end
        if ~isfield(params,'TauStd'); params.TauStd = defparams.TauStd; end
        if ~isfield(params,'A_lb'); params.A_lb = defparams.A_lb; end
        if ~isfield(params,'b_lb'); params.b_lb = defparams.b_lb; end
        if ~isfield(params,'c1_lb'); params.c1_lb = defparams.c1_lb; end
        if ~isfield(params,'print_flag'); params.print_flag = defparams.print_flag; end
    end
    
    %% ZW: compute SAM (no MCMC result as initialization)
    Dt = 1;                                     % length of time bin
    marg_flag = params.marg;
    gam_flag = params.upd_gam;
    gam_step = params.gam_step;
    std_move = params.std_move;
    add_move = params.add_move;
    if isempty(params.g)
        p = params.p;
    else
        p = length(params.g);                       % order of autoregressive process
    end
    if isempty(params.init)
       fprintf('Initializing using noise constrained FOOPSI...  ');
       params.init = get_initial_sample(Y,params);
       fprintf('done. \n');
    end 
    SAM = params.init;
    
    %% ZW: Try if MCMC works
    try
        g = SAM.g(:)';
        if g == 0 %#ok<BDSCI>
            gr = [0.9,0.1];
            pl = poly(gr);
            g = -pl(2:end);
            p = 2;
        end
        gr = sort(roots([1,-g(:)']));
        if p == 1; gr = [0,gr]; end
        if any(gr<0) || any(~isreal(gr)) || length(gr)>2 || max(gr)>0.998
            gr = params.defg;
        end
        tau = -Dt./log(gr);
        tau1_std = max(tau(1)/100,params.TauStd(1));
        tau2_std = min(tau(2)/5,params.TauStd(2)); 
        ge = max(gr).^(0:T-1)';
        if p == 1
            G1 = sparse(1:T,1:T,Inf*ones(T,1));
        elseif p == 2
            G1 = spdiags(ones(T,1)*[-min(gr),1], -1:0,T,T);
        else
            error('This order of the AR process is currently not supported');
        end
        G2 = spdiags(ones(T,1)*[-max(gr),1], -1:0,T,T);
        sg = SAM.sg;
        SAM = params.init;    
        spiketimes_ = SAM.spiketimes_;
        lam_ = SAM.lam_;
        A_ = SAM.A_*diff(gr);
        b_ = SAM.b_;
        C_in = SAM.C_in;        
        s_1 = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau(1)),T,1);  
        s_2 = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau(2)),T,1);  
        if ~isfield(params,'prec') 
            % FN --
            % Eftychios: prec specifies to what extent you want to discard 
            % the long slowly decaying tales of the ca response. Try setting 
            % it e.g., to 5e-2 instead of 1e-2 to speed things up.
            prec = 1e-2;     % precision
        else
            prec = params.prec;
        end
        ef_d = exp(-(0:T)/tau(2));
        if p == 1
            h_max = 1; % max value of transient    
            ef_h = [0,0];
            e_support = find(ef_d<prec*h_max,1,'first');
            if isempty(e_support);
                e_support = T;
            end
            e_support = min(e_support,T);
        else
            t_max = (tau(1)*tau(2))/(tau(2)-tau(1))*log(tau(2)/tau(1)); %time of maximum
            h_max = exp(-t_max/tau(2)) - exp(-t_max/tau(1)); % max value of transient    
            ef_h = -exp(-(0:T)/tau(1));
            e_support = find(ef_d-ef_h<prec*h_max,1,'first');
            if isempty(e_support);
                e_support = T;
            end
            e_support = min(e_support,T);
        end
        ef_h = ef_h(1:min(e_support,length(ef_h)))/diff(gr);
        ef_d = ef_d(1:e_support)/diff(gr);
        ef = [{ef_h ef_d};{cumsum(ef_h.^2) cumsum(ef_d.^2)}];
        B = params.B;
        N = params.Nsamples + B;
        if p == 1; G1sp = zeros(T,1); else G1sp = G1\s_1(:); end
        Gs = (-G1sp(:)+G2\s_2(:))/diff(gr);
        ss = cell(N,1); 
        lam = zeros(N,1);
        Am = zeros(N,1);
        ns = zeros(N,1);
        Gam = zeros(N,2);
        if ~marg_flag
            Cb = zeros(N,1);
            Cin = zeros(N,1);
            SG = zeros(N,1);
        end
        Sp = .1*range(Y)*eye(3);  % prior covariance
        Ld = inv(Sp);
        lb = [params.A_lb/h_max*diff(gr),params.b_lb,params.c1_lb]'; % lower bound for [A,Cb,Cin]
        A_ = max(A_,1.1*lb(1));
        mu = [A_;b_;C_in]; % prior mean 
        Ns = 15; % Number of HMC samples
        Ym = Y - ones(T,1)*mu(2) - ge*mu(3);
        mub = zeros(2,1);
        Sigb = zeros(2,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extra tau-related params
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        tauMoves = [0 0];
        tau_min = 0;
        tau_max = 500;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:N
            if gam_flag
                Gam(i,:) = tau;
            end
            sg_ = sg;
            rate = @(t) lambda_rate(t,lam_);
            [spiketimes, ~]  = get_next_spikes(spiketimes_(:)',A_*Gs',Ym',ef,tau,sg_^2, rate, std_move, add_move, Dt, A_);
            spiketimes_ = spiketimes;
            spiketimes(spiketimes<0) = -spiketimes(spiketimes<0);
            spiketimes(spiketimes>T*Dt) = 2*T*Dt - spiketimes(spiketimes>T*Dt); 
            trunc_spikes = ceil(spiketimes/Dt);
            trunc_spikes(trunc_spikes == 0) = 1;
            s_1 =   sparse(trunc_spikes,1,exp((spiketimes_ - Dt*trunc_spikes)/tau(1)),T,1);
            s_2 =   sparse(trunc_spikes,1,exp((spiketimes_ - Dt*trunc_spikes)/tau(2)),T,1);  
            if p == 1; G1sp = zeros(T,1); else G1sp = G1\s_1(:); end
            Gs = (-G1sp+G2\s_2(:))/diff(gr);
            ss{i} = spiketimes;
            nsp = length(spiketimes);
            ns(i) = nsp;
            lam(i) = nsp/(T*Dt);
            lam_ = lam(i);    
            AM = [Gs,ones(T,1),ge];
            L = inv(Ld + AM'*AM/sg^2);
            mu_post = (Ld + AM'*AM/sg^2)\(AM'*Y/sg^2 + Sp\mu);
            if ~marg_flag
                x_in = [A_;b_;C_in];
                if any(x_in < lb)
                    x_in = max(x_in,1.1*lb);
                end
                if all(isnan(L(:))) 
                    % FN --
                    % added to avoid error in R = chol(L) in HMC_exact2 due 
                    % to L not being positive definite. It happens when 
                    % isnan(det(Ld + AM'*AM/sg^2)), ie when Ld + AM'*AM/sg^2 
                    % is singular (not invertible).
                    Am(i) = NaN;
                    Cb(i) = NaN;
                    Cin(i) = NaN';
                else        
                    [temp,~] = HMC_exact2(eye(3), -lb, L, mu_post, 1, Ns, x_in);
                    Am(i) = temp(1,Ns);
                    Cb(i) = temp(2,Ns);
                    Cin(i) = temp(3,Ns)';
                end
                A_ = Am(i);
                b_ = Cb(i);
                C_in = Cin(i);

                Ym   = Y - b_ - ge*C_in;
                res   = Ym - A_*Gs;
                sg   = 1./sqrt(gamrnd(1+T/2,1/(0.1 + sum((res.^2)/2))));
                SG(i) = sg;
            else
                repeat = 1;
                while repeat
                    A_ = mu_post(1) + sqrt(L(1,1))*randn;
                    repeat = (A_<0);
                end                
                Am(i) = A_;
                if i > B
                   mub = mub + mu_post(2+(0:p));
                   Sigb = Sigb + L(2+(0:p),2+(0:p));
                end
            end
            if gam_flag
                if mod(i-B,gam_step) == 0  % update time constants
                    if p >= 2       % update rise time constant
                        logC = -norm(Y(:) - AM*[A_;b_;C_in])^2; 
                        tau_ = tau;
                        tau_temp = tau_(1)+(tau1_std*randn); 
                        while tau_temp >tau(2) || tau_temp<tau_min
                            tau_temp = tau_(1)+(tau1_std*randn);
                        end 
                        tau_(1) = tau_temp;
                        gr_ = exp(Dt*(-1./tau_));
                        s_1_ =   sparse(trunc_spikes,1,exp((spiketimes_ - Dt*trunc_spikes)/tau_(1)),T,1);  
                        G1_ = spdiags(ones(T,1)*[-min(gr_),1],[-1:0],T,T);
                        Gs_ = (-G1_\s_1_(:)+G2\s_2(:))/diff(gr_);

                        logC_ = -norm(Y(:)-A_*Gs-b_-C_in*ge)^2;
                        %accept or reject
                        prior_ratio = 1;
                %        prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
                        ratio = exp((logC_-logC)/(2*sg^2))*prior_ratio;
                        if rand < ratio %accept
                            tau = tau_;
                            G1 = G1_; %c = c_; 
                            Gs = Gs_; gr = gr_; s_1 = s_1_;
                            tauMoves = tauMoves + [1 1];
                        else
                            tauMoves = tauMoves + [0 1];
                        end                
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%
                    % next update decay time constant
                    %%%%%%%%%%%%%%%%%%%%%%%
                    %initial logC
                    logC = -norm(Y(:)-A_*Gs-b_-C_in*ge)^2;
                    tau_ = tau;
                    tau_temp = tau_(2)+(tau2_std*randn);
                    while tau_temp>tau_max || tau_temp<tau_(1)
                        tau_temp = tau_(2)+(tau2_std*randn);
                    end  
                    tau_(2) = tau_temp;
                    s_2_ =   sparse(trunc_spikes,1,exp((spiketimes_ - Dt*trunc_spikes)/tau_(2)),T,1);  
                    gr_ = exp(Dt*(-1./tau_));
                    ge_ = max(gr_).^(0:T-1)';
                    G2_ = spdiags(ones(T,1)*[-max(gr_),1],[-1:0],T,T);  
                    if p == 1; G1sp = zeros(T,1); else G1sp = G1\s_1(:); end
                    Gs_ = (-G1sp+G2_\s_2_(:))/diff(gr_);
                    logC_ = -norm(Y(:)-A_*Gs_-b_-C_in*ge_)^2;
                    %accept or reject
                    prior_ratio = 1;
            %       prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
                    ratio = exp((1./(2*sg^2)).*(logC_-logC))*prior_ratio;
                    if rand<ratio %accept
                        tau = tau_;
                        %c = c_; 
                        ge = ge_; Gs = Gs_; G2 = G2_; gr = gr_; s_2 = s_2_; %#ok<NASGU>
                        tauMoves = tauMoves + [1 1];
                    else
                        tauMoves = tauMoves + [0 1];
                    end                
                    ef_d = exp(-(0:T)/tau(2));
                    if p == 1
                        h_max = 1; % max value of transient    
                        ef_h = [0,0];
                        e_support = find(ef_d<prec*h_max,1,'first');
                        if isempty(e_support);
                            e_support = T;
                        end
                        e_support = min(e_support,T);
                    else
                        t_max = (tau(1)*tau(2))/(tau(2)-tau(1))*log(tau(2)/tau(1)); %time of maximum
                        h_max = exp(-t_max/tau(2)) - exp(-t_max/tau(1)); % max value of transient    
                        ef_h = -exp(-(0:T)/tau(1));
                        e_support = find(ef_d-ef_h<prec*h_max,1,'first');
                        if isempty(e_support);
                            e_support = T;
                        end
                        e_support = min(e_support,T);
                    end
                    ef_h = ef_h(1:min(e_support,length(ef_h)))/diff(gr);
                    ef_d = ef_d(1:e_support)/diff(gr);
                    ef = [{ef_h ef_d};{cumsum(ef_h.^2) cumsum(ef_d.^2)}];     
                end
            end
            if params.print_flag && mod(i,100)==0
                fprintf('%i out of total %i samples drawn \n', i, N);
            end 
        end
        if marg_flag
            mub = mub/(N-B);
            Sigb = Sigb/(N-B)^2;
        end

        if marg_flag
            SAMPLES.Cb = [mub(1),sqrt(Sigb(1,1))];
            SAMPLES.Cin = [mub(1+(1:p)),sqrt(diag(Sigb(1+(1:p),1+(1:p))))];
        else
            SAMPLES.Cb = Cb(B+1:N);
            SAMPLES.Cin = Cin(B+1:N,:);
            SAMPLES.sn2 = SG(B+1:N).^2;
        end
        SAMPLES.ns = ns(B+1:N);
        SAMPLES.ss = ss(B+1:N);
        SAMPLES.ld = lam(B+1:N);
        SAMPLES.Am = Am(B+1:N);
        if gam_flag
            SAMPLES.g = Gam(B+1:N,:);
        end
        SAMPLES.params = params.init; % this is SAM
        SAMPLES.C_rec  = make_mean_sample(SAMPLES,Y);
        SAMPLES.F_est  = mean(SAMPLES.C_rec);
        SAMPLES.spikeRaster = samples_cell2mat(SAMPLES.ss,T);
        SAMPLES.spk    = mean(SAMPLES.spikeRaster);
        SAMPLES.HMC2   = true;
    catch
        SAMPLES        = SAM;
        SAMPLES.params = params.init;
        SAMPLES.HMC2   = false;
    end
end

function SAM = get_initial_sample(Y,params)
    % obtain initial sample by performing sparse noise-constrained deconvolution
    % Author: Eftychios A. Pnevmatikakis    
    if isfield(params,'p'); options.p = params.p; else options.p = 1; end
    [c,b,c1,g,sn,sp] = constrained_foopsi(Y,params.b,params.c1,params.g,params.sn,options);
    Dt = 1;
    T = length(Y);
    if ~exist('sp','var')
        G = make_G_matrix(T,params.g);
        sp = G*c;
    end
    s_in = sp>0.15*max(sp);
    spiketimes_ = Dt*(find(s_in) + rand(size(find(s_in))) - 0.5);
    spiketimes_(spiketimes_ >= T*Dt) = 2*T*Dt - spiketimes_(spiketimes_ >= T*Dt);
    SAM.c    = c;
    SAM.F_est = c;
    SAM.sp   = sp;
    SAM.spk  = sp;
    SAM.spk(~s_in) = 0;
    SAM.s_in = s_in;
    SAM.lam_ = length(spiketimes_)/(T*Dt);
    SAM.spiketimes_ = spiketimes_;
    SAM.A_   = max(median(sp(s_in)),max(sp(s_in))/4);  % initial amplitude value
    if length(g) == 2
        SAM.A_   = SAM.A_/sqrt(g(1)^2+4*g(2));
    end
    SAM.b_   = max(b,range(Y)/25);                     % initial baseline value
    SAM.C_in = max(c1,(Y(1)-b)/10);                    % initial value sample
    SAM.sg = sn;                                       % initial noise value
    SAM.g = g;                                         % initial time constant value
end

function [samples, ci]  = get_next_spikes(curr_spikes,curr_calcium,calciumSignal,ef,tau,calciumNoiseVar, lam, proposalVar, add_move, Dt, A)
    %addMoves, dropMoves, and timeMoves give acceptance probabilities for each subclass of move
    %the samples will be a cell array of lists of spike times - the spike times won't be sorted but this shouldn't be a problem.
    %noise level here matters for the proposal distribution (how much it 
    %should trust data for proposals vs how much it should mix based on uniform prior)
    %this is accounted for by calciumNoiseVar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initialize some parameters
    T = length(calciumSignal); %for all of this, units are bins and spiketrains go from 0 to T where T is number of bins
    %   nsweeps = 1e3; %number of sweeps.
    nsweeps = 1;
    samples = cell(nsweeps);
    %% start with initial spiketrain and initial predicted calcium 
    si = curr_spikes; %initial set of spike times has no spikes - this will not be sorted but that shouldn't be a problem
    ci = curr_calcium; %initial calcium is set to baseline 
    
    N = length(si); %number of spikes in spiketrain
        
    %initial logC - compute likelihood initially completely - updates to likelihood will be local
    logC = -(ci-calciumSignal)*(ci-calciumSignal)'; 
    %m = p_spike*Dt*T;    
    %flag for uniform vs likelihood proposal (if using likelihood proposal, then time shifts are pure Gibbs)
    %this really should be split into four cases, 
    % 1) RW for time shifts with uniform add/drop
    % 2) RW for time shifts with likelihood proposal add/drop
    % 3) Gibbs for time shifts with uniform add/drop
    % 4) Gibbs for time shifts with likelihood proposal add/drop    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% loop over sweeps to generate samples
    addMoves = [0 0]; %first elem is number successful, second is number total
    dropMoves = [0 0];
    timeMoves = [0 0];
    time_move = 0;
    time_add = 0;
    for i = 1:nsweeps
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% loop over spikes, perform spike time move (could parallelize here for non-interacting spikes, i.e. spikes that are far enough away)
        %move
        for ni = 1:N %possibly go through in a random order (if you like)
            tmpi = si(ni);
            tmpi_ = si(ni)+(proposalVar*randn); %with bouncing off edges
            if tmpi_<0
                tmpi_ = -(tmpi_);
            elseif tmpi_>T
                tmpi_ = T-(tmpi_-T);
            end
            %set si_ to set of spikes with the move and ci_ to adjusted calcium and update logC_ to adjusted
            [si_, ci_, logC_] = removeSpike(si,ci,logC,ef,tau,calciumSignal,tmpi,ni,Dt,A);
            [si_, ci_, logC_] = addSpike(si_,ci_,logC_,ef,tau,calciumSignal,tmpi_,Dt,A);
            %accept or reject
            ratio = exp((1/(2*calciumNoiseVar))*(logC_-logC)*lam(tmpi)/lam(tmpi_));
            if ratio>1 %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                timeMoves = timeMoves + [1 1];
            elseif rand<ratio %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                timeMoves = timeMoves + [1 1];
            else
                %reject - do nothing
                timeMoves = timeMoves + [0 1];
            end
        end
        N = length(si);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% loop over add/drop a few times
        %define insertion proposal distribution as the likelihood function
        %define removal proposal distribution as uniform over spikes
        %perhaps better is to choose smarter removals.
        for ii = 1:add_move 
            %% add
            %propose a uniform add
            tmpi = T*Dt*rand;         
            [si_, ci_, logC_] = addSpike(si,ci,logC,ef,tau,calciumSignal,tmpi,Dt,A);
            %forward probability
            fprob = 1/(T*Dt);
            %reverse (remove at that spot) probability
            rprob = 1/(N+1);
            %accept or reject
            ratio = exp((1/(2*calciumNoiseVar))*(logC_ - logC))*(rprob/fprob)*lam(tmpi); %posterior times reverse prob/forward prob
            if ratio>1 %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                addMoves = addMoves + [1 1];
            elseif rand<ratio %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                addMoves = addMoves + [1 1];
            else
                %reject - do nothing
                addMoves = addMoves + [0 1];
            end
            N = length(si);

            %% delete
            if N>0                
                %propose a uniform removal
                tmpi = randi(N);
                [si_, ci_, logC_] = removeSpike(si,ci,logC,ef,tau,calciumSignal,si(tmpi),tmpi,Dt,A);
                %reverse probability
                rprob = 1/(T*Dt);
                %compute forward prob
                fprob = 1/N;                
                %accept or reject
                ratio = exp((1/(2*calciumNoiseVar))*(logC_ - logC))*(rprob/fprob)*(1/lam(si(tmpi))); %posterior times reverse prob/forward prob
                if ratio>1 %accept
                    si = si_;
                    ci = ci_;
                    logC = logC_;
                    dropMoves = dropMoves + [1 1];
                elseif rand<ratio %accept
                    si = si_;
                    ci = ci_;
                    logC = logC_;
                    dropMoves = dropMoves + [1 1];
                else
                    %reject - do nothing
                    dropMoves = dropMoves + [0 1];
                end
                N = length(si);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        N_sto = [N_sto N];
        samples = si;
        %store overall logliklihood as well
    %        objective = [objective logC];
    end
end

function [Xs, bounce_count] = HMC_exact2(F, g, M, mu_r, cov, L, initial_X)
    % Author: Ari Pakman
    % returns samples from a d-dimensional Gaussian with m constraints given by  F*X+g >0 
    % If cov == true
    % then M is the covariance and the mean is mu = mu_r 
    % if cov== false 
    % then M is the precision matrix and the log-density is -1/2 X'*M*X + r'*X
    % Input
    % F:          m x d matrix
    % g:          m x 1 vector 
    % M           d x d matrix, must be symmmetric and definite positive
    % mu_r        d x 1 vector. 
    % cov:        see explanation above 
    % L:          number of samples desired
    % initial_X   d x 1 vector. Must satisfy the constraint.
    % Output
    % Xs:      d x L matrix, each column is a sample
    % bounce_count:  number of times the particle bounced 
    %% go to a whitened frame
    m = size(g,1);
    if size(F,1) ~= m
        display('error');
        return
    end
    if cov 
        mu= mu_r;
        g = g + F*mu;
        %if min(eig(M))<0
        %    M = M - 1.01*min(eig(M))*eye(size(M,1));
        %end
        R = chol(M);
        F = F*R';
        initial_X= initial_X -mu;
        initial_X = R'\initial_X;
    else
        r=mu_r;
        R=chol(M);      % M = R'*R    
        mu = R\(R'\r);
        g = g+F*mu;
        F = F/R;               % this is the time-consuming step in this code section
        initial_X= initial_X -mu;
        initial_X = R*initial_X;    
    end
    d = size(initial_X,1);
    Xs=NaN;
    bounce_count =0;
    nearzero= 10000*eps;
    % Verify that initial_X is feasible
     c= F*initial_X +g;
     if any(c<0)
         display('error: inconsistent initial condition');
         return
     end
    % Unsparcify the linear constraints
    g=full(g);
    F2 = sum(F.*F,2);  % squared norm of the rows of F, needed for reflecting the velocity
    F=full(F);         % if we don't unsparcify  qj= F(j,:)*V/F2(j) becomes very slow.
    Ft = F';
    % Sampling loop
    last_X= initial_X;
    Xs=zeros(d,L);
    Xs(:,1)=initial_X;
    i=2;
    while (i <= L)
    %i
    stop=0;   
    j=0;
    V0= normrnd(0,1, d,1);   % initial velocity
    X = last_X;
    T=pi/2;                  % total time the particle will move
    tt=0;                    % records how much time the particle already moved 
        while (1) 
            a = V0; 
            a= real(a);
            b = X;
            fa = F*a;
            fb = F*b;
            U = sqrt(fa.^2 + fb.^2);
            phi = atan2(-fa,fb);           % -pi < phi < +pi    
            pn = abs(g./U)<=1; % these are the walls that may be hit 
            % find the first time constraint becomes zero
            if any(pn) 
                inds = find(pn);
                phn= phi(pn);
                t1=-phn + acos(-g(pn)./U(pn));  % time at which coordinates hit the walls                                                 % t1 in [-pi, 2*pi]
                t1(t1<0) = 2*pi + t1(t1<0);     % t1 in [0, 2*pi]                          
                t2 = -t1 -2*phn;                % second solution to hit the walls 
                t2(t2<0) = 2*pi + t2(t2<0);     % t2 in [-2*pi, 2*pi] 
                t2(t2<0) = 2*pi + t2(t2<0);     % t2 in [0, 2*pi] 
             % if there was a previous reflection (j>0)
             % and there is a potential reflection at the sample plane                                    
             % make sure that a new reflection at j is not found because of numerical error
                if j>0    
                    if pn(j)==1    
                        indj=sum(pn(1:j));            
                        tt1 = t1(indj);
                        if abs(tt1) < nearzero || abs(tt1-2*pi)< nearzero
                            t1(indj)=Inf;
                        else
                            tt2 = t2(indj);
                            if abs(tt2) < nearzero || abs(tt1-2*pi)< nearzero 
                                t2(indj) = Inf;
                            end
                        end                    
                    end
                end
                [mt1, ind1] = min(t1);
                [mt2, ind2] = min(t2);
                [mt, ind12]  = min([mt1 mt2]);
                if ind12==1
                    m_ind = ind1;
                else
                    m_ind= ind2;
                end
                % find the reflection plane 
                j = inds(m_ind);      % j is an index in the full vector of dim-m, not in the restriced vector determined by pn.
            else  %if pn(i) =0 for all i
                    mt =T;
            end   %if pn(i)
            tt=tt+mt;
            if tt>=T
                mt= mt-(tt-T);
                stop=1;
            end
            % move the particle a time mt
            X = a*sin(mt) + b*cos(mt);
            V = a*cos(mt) - b*sin(mt);
            if stop                    
                break;
            end
            % compute reflected velocity
            qj=  F(j,:)*V/F2(j);   
            V0 = V -2*qj*Ft(:,j);
            bounce_count = bounce_count+ 1;
            % dif =V0'*M*V0 - V'*M*V;
        end % while(1)
        % at this point we have a sampled value X, but due to possible
        % numerical instabilities we check that the candidate X satisfies the
        % constraints before accepting it.
        if all(F*X +g > 0)
            Xs(:,i)=X;
            last_X = X;
            i= i+1;
        else
            % disp('hmc reject')
            % ZW -- 
            error('hmc reject')
        end 
    end %while (i <= L)
    % transform back to the unwhitened frame
    if cov
        Xs = R'*Xs  + repmat(mu, 1,L);   
    else
        Xs = R\Xs  + repmat(mu, 1,L);   
    end
end

function lam = lambda_rate(t,lambda) %#ok<INUSL>
    lam = lambda;
end

function G = make_G_matrix(T,g,varargin)
    if length(g) == 1 && g < 0
        g=0;
    end
    G = spdiags(ones(T,1)*[-flipud(g(:))',1],-length(g):0,T,T);
    if nargin == 3
        sl = [0;cumsum(varargin{1}(:))];
        for i = 1:length(sl)-1
            G(sl(i)+1,sl(i+1))=0;
        end
    end
end

function [newSpikeTrain, newCalcium, newLL] = removeSpike(oldSpikeTrain,oldCalcium,oldLL,filter,tau,obsCalcium,timeToRemove,indx,Dt,A)
    tau_h = tau(1);
    tau_d = tau(2);
    ef_h = filter{1,1};
    ef_d = filter{1,2};
    ef_nh = filter{2,1};
    ef_nd = filter{2,2};
    newSpikeTrain = oldSpikeTrain;
    newSpikeTrain(indx) = [];
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient    
    wk_h = A*exp((timeToRemove - Dt*ceil(timeToRemove/Dt))/tau_h);
    wk_d = A*exp((timeToRemove - Dt*ceil(timeToRemove/Dt))/tau_d);
    %handle ef_h first
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(ef_h)+floor(timeToRemove)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) - wk_h*ef_h(1:length(tmp));
    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    %newLL = oldLL - ( wk_h^2*norm(ef_h(1:length(tmp)))^2 + 2*relevantResidual*(wk_h*ef_h(1:length(tmp))'));
    newLL = oldLL - ( wk_h^2*ef_nh(length(tmp)) + 2*relevantResidual*(wk_h*ef_h(1:length(tmp))'));
    oldCalcium = newCalcium;
    oldLL = newLL;
    %handle ef_d next
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(ef_d)+floor(timeToRemove)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) - wk_d*ef_d(1:length(tmp));
    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    %newLL = oldLL - ( wk_d^2*norm(ef_d(1:length(tmp)))^2 + 2*relevantResidual*(wk_d*ef_d(1:length(tmp))'));
    newLL = oldLL - ( wk_d^2*ef_nd(length(tmp)) + 2*relevantResidual*(wk_d*ef_d(1:length(tmp))'));
end

function [newSpikeTrain, newCalcium, newLL] = addSpike(oldSpikeTrain,oldCalcium,oldLL,filter,tau,obsCalcium,timeToAdd,Dt,A)
    tau_h = tau(1);
    tau_d = tau(2);
    ef_h = filter{1,1};
    ef_d = filter{1,2};
    ef_nh = filter{2,1};
    ef_nd = filter{2,2};
    newSpikeTrain = [oldSpikeTrain timeToAdd]; %possibly inefficient, change if problematic (only likely to be a problem for large numbers of spikes)
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient
    wk_h = A*exp((timeToAdd - Dt*ceil(timeToAdd/Dt))/tau_h);
    wk_d = A*exp((timeToAdd - Dt*ceil(timeToAdd/Dt))/tau_d);    
    %%%%%%%%%%%%%%%%%
    %handle ef_h first
    newCalcium = oldCalcium;
    tmp = 1 + (floor(timeToAdd):min((length(ef_h)+floor(timeToAdd)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) + wk_h*ef_h(1:length(tmp));
    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    %newLL = oldLL - ( wk_h^2*norm(ef_h(1:length(tmp)))^2 - 2*relevantResidual*(wk_h*ef_h(1:length(tmp))'));
    newLL = oldLL - ( wk_h^2*ef_nh(length(tmp)) - 2*relevantResidual*(wk_h*ef_h(1:length(tmp))'));
    oldCalcium = newCalcium;
    oldLL = newLL;
    %handle ef_d next
    tmp = 1 + (floor(timeToAdd):min((length(ef_d)+floor(timeToAdd)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) + wk_d*ef_d(1:length(tmp));
    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    %newLL = oldLL - ( wk_d^2*norm(ef_d(1:length(tmp)))^2 - 2*relevantResidual*(wk_d*ef_d(1:length(tmp))'));
    newLL = oldLL - ( wk_d^2*ef_nd(length(tmp)) - 2*relevantResidual*(wk_d*ef_d(1:length(tmp))'));
end


function C_rec = make_mean_sample(SAMPLES,Y)
    T = length(Y);
    N = length(SAMPLES.ns);
    show_gamma = 0; %#ok<NASGU>
    P = SAMPLES.params;
    P.f = 1;
    g = P.g(:);
    p = length(g); %#ok<NASGU>
    Dt = 1/P.f;                                     % length of time bin
    if ~isfield(SAMPLES,'g');
        show_gamma = 0; %#ok<NASGU>
        SAMPLES.g = ones(N,1)*g';
    end
    if length(SAMPLES.Cb) == 2
        marg = 1;       % marginalized sampler
    else
        marg = 0;       % full sampler
    end
    C_rec = zeros(N,T);
    for rep = 1:N
        %trunc_spikes = ceil(SAMPLES.ss{rep}/Dt);
        tau = SAMPLES.g(rep,:);
        gr = exp(-1./tau);
        ge = max(gr).^(0:T-1)';
        s_1 =   sparse(ceil(SAMPLES.ss{rep}/Dt),1,exp((SAMPLES.ss{rep} - Dt*ceil(SAMPLES.ss{rep}/Dt))/tau(1)),T,1);  
        s_2 =   sparse(ceil(SAMPLES.ss{rep}/Dt),1,exp((SAMPLES.ss{rep} - Dt*ceil(SAMPLES.ss{rep}/Dt))/tau(2)),T,1);  
        if gr(1) == 0
            G1 = sparse(1:T,1:T,Inf*ones(T,1)); %#ok<NASGU>
            G1sp = zeros(T,1);
        else
            G1 = spdiags(ones(T,1)*[-min(gr),1],[-1:0],T,T);
            G1sp = G1\s_1(:);
        end    
        G2 = spdiags(ones(T,1)*[-max(gr),1],[-1:0],T,T);
        Gs = (-G1sp+ G2\s_2(:))/diff(gr);
        if marg
            %C_rec(rep,:) = SAMPLES.Cb(1) + SAMPLES.Am(rep)*filter(1,[1,-SAMPLES.g(rep,:)],full(s_)+[SAMPLES.Cin(:,1)',zeros(1,T-p)]);
            C_rec(rep,:) = SAMPLES.Cb(1) + SAMPLES.Am(rep)*Gs + (ge*SAMPLES.Cin(:,1));
        else
            %C_rec(rep,:) = SAMPLES.Cb(rep) + SAMPLES.Am(rep)*filter(1,[1,-SAMPLES.g(rep,:)],full(s_)+[SAMPLES.Cin(rep,:),zeros(1,T-p)]);
            C_rec(rep,:) = SAMPLES.Cb(rep) + SAMPLES.Am(rep)*Gs + (ge*SAMPLES.Cin(rep,:)');
        end
    end
end

function spikeRaster = samples_cell2mat(sampleCell,T,Dt)
    if nargin == 2
        Dt = 1;
    end
    bins = 0:Dt:(T-Dt);
    nsamples = length(sampleCell);
    spikeRaster = zeros(nsamples,length(bins));
    for i = 1:nsamples
        tmp = histc([sampleCell{i}(:); inf],[bins, (T+1)]);
        spikeRaster(i,:) = tmp(1:(end-1))';
    end
end


% function [samples, addMoves, dropMoves, timeMoves, N_sto]  = sampleSpikes(calciumSignal,ef,tau,b,calciumNoiseVar, p_spike, proposalVar, nsweeps)
%     % ZW --
%     % This function is somewhat identical to get_next_spike
%     %
%     %
%     %
%     %addMoves, dropMoves, and timeMoves give acceptance probabilities for each subclass of move
%     %the samples will be a cell array of lists of spike times - the spike times won't be sorted but this shouldn't be a problem.
%     %noise level here matters for the proposal distribution (how much it 
%     %should trust data for proposals vs how much it should mix based on uniform prior)
%     %this is accounted for by calciumNoiseVar    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % initialize some parameters
%     T = length(calciumSignal); %for all of this, units are bins and spiketrains go from 0 to T where T is number of bins
%     % nsweeps = 1e3; %number of sweeps.
%     samples = cell(nsweeps);
%     N_sto = [];
%     objective = [];
%     % start with initial spiketrain and initial predicted calcium 
%     si = []; %initial set of spike times has no spikes - this will not be sorted but that shouldn't be a problem
%     ci = b*ones(1,T); %initial calcium is set to baseline 
%     N = length(si); %number of spikes in spiketrain
%     %initial logC - compute likelihood initially completely - updates to likelihood will be local
%     logC = -(ci-calciumSignal)*(ci-calciumSignal)'; 
%     m = p_spike*T;
%     %flag for uniform vs likelihood proposal (if using likelihood proposal, then time shifts are pure Gibbs)
%     %this really should be split into four cases, 
%     % 1) RW for time shifts with uniform add/drop
%     % 2) RW for time shifts with likelihood proposal add/drop
%     % 3) Gibbs for time shifts with uniform add/drop
%     % 4) Gibbs for time shifts with likelihood proposal add/drop
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% loop over sweeps to generate samples
%     addMoves = [0 0]; %first elem is number successful, second is number total
%     dropMoves = [0 0];
%     timeMoves = [0 0];
%     for i = 1:nsweeps
%         % loop over spikes, perform spike time move (could parallelize here for non-interacting spikes, i.e. spikes that are far enough away)
%         %move
%         for ni = 1:N %possibly go through in a random order (if you like)
%             tmpi = si(ni);
%             tmpi_ = si(ni)+(proposalVar*randn); %with bouncing off edges
%             if tmpi_<0
%                 tmpi_ = -(tmpi_);
%             elseif tmpi_>T
%                 tmpi_ = T-(tmpi_-T);
%             end
%             %set si_ to set of spikes with the move and ci_ to adjusted calcium and update logC_ to adjusted
%             [si_, ci_, logC_] = removeSpike(si,ci,logC,ef,tau,calciumSignal,tmpi,ni);
%             [si_, ci_, logC_] = addSpike(si_,ci_,logC_,ef,tau,calciumSignal,tmpi_);
%             %accept or reject
%             ratio = exp((1/(2*calciumNoiseVar))*(logC_-logC));
%             if ratio>1 %accept
%                 si = si_;
%                 ci = ci_;
%                 logC = logC_;
%                 timeMoves = timeMoves + [1 1];
%             elseif rand<ratio %accept
%                 si = si_;
%                 ci = ci_;
%                 logC = logC_;
%                 timeMoves = timeMoves + [1 1];
%             else
%                 %reject - do nothing
%                 timeMoves = timeMoves + [0 1];
%             end
%         end        
%         N = length(si);
%         % loop over add/drop a few times
%         %define insertion proposal distribution as the likelihood function
%         %define removal proposal distribution as uniform over spikes
%         %perhaps better is to choose smarter removals.
%         for ii = 1:10 
%             % add
%             %propose a uniform add
%             tmpi = T*rand;         
%             [si_, ci_, logC_] = addSpike(si,ci,logC,ef,tau,calciumSignal,tmpi);
%             %forward probability
%             fprob = 1/T;
%             %reverse (remove at that spot) probability
%             rprob = 1/(N+1);
%             %accept or reject
%             ratio = exp((1/(2*calciumNoiseVar))*(logC_ - logC))*(rprob/fprob)*(m/(T-m)); %posterior times reverse prob/forward prob
%             if ratio>1 %accept
%                 si = si_;
%                 ci = ci_;
%                 logC = logC_;
%                 addMoves = addMoves + [1 1];
%             elseif rand<ratio %accept
%                 si = si_;
%                 ci = ci_;
%                 logC = logC_;
%                 addMoves = addMoves + [1 1];
%             else
%                 %reject - do nothing
%                 addMoves = addMoves + [0 1];
%             end
%             N = length(si);
%             % delete
%             if N>0                
%                 %propose a uniform removal
%                 tmpi = randi(N);
%                 [si_, ci_, logC_] = removeSpike(si,ci,logC,ef,tau,calciumSignal,si(tmpi),tmpi);
%                 %reverse probability
%                 rprob = 1/T;
%                 %compute forward prob
%                 fprob = 1/N;
%                 %accept or reject
%                 ratio = exp((1/(2*calciumNoiseVar))*(logC_ - logC))*(rprob/fprob)*((T-m)/m); %posterior times reverse prob/forward prob
%                 if ratio>1 %accept
%                     si = si_;
%                     ci = ci_;
%                     logC = logC_;
%                     dropMoves = dropMoves + [1 1];
%                 elseif rand<ratio %accept
%                     si = si_;
%                     ci = ci_;
%                     logC = logC_;
%                     dropMoves = dropMoves + [1 1];
%                 else
%                     %reject - do nothing
%                     dropMoves = dropMoves + [0 1];
%                 end
%                 N = length(si);
%             
%             end
%         end
%         N_sto = [N_sto N];
%         samples{i} = si;
%         %store overall logliklihood as well
%         objective = [objective logC];        
%         % figure(48)
%         % subplot(211)
%         % plot(N_sto)
%         % subplot(212)
%         % plot(objective)
%         
%         % [addMoves(1)/addMoves(2) dropMoves(1)/dropMoves(2)]
%     end
% end

% function [f,grad] = min_gamma(g,s,y)
%     T = length(s);
%     prec = 1e-4;
%     if length(g) == 1
%         K = min(ceil(log(prec)/log(g)),T);
%         vec_f = g.^(0:K-1);
%         G_f = toeplitz(sparse(1:K,1,vec_f(:),T,1),sparse(1,1,1,1,T));
%         vec_g = (0:K-1).*(g.^(-1:K-2));
%         G_g = toeplitz(sparse(1:K,1,vec_g(:),T,1),sparse(1,T));
% 
%         f = 0.5*norm(y-G_f*s)^2;
%         grad = -(y-G_f*s)'*(G_g*s);
%     else
%         G = spdiags(ones(T,1)*[-flipud(g(:))',1],-length(g):0,T,T);
%         Gi = G\[1;zeros(T-1,1)];
%         K = min(find(Gi<prec,1,'first'),T);
%         Gs = G\s;
%         f = 0.5*norm(y-Gs)^2;
%         g1 = g(1); g2 = g(2);
%         d = sqrt(g1^2+4*g2);
%         v1 = 2.^(1:K).*(g1*((g1-d).^(1:K) - (g1+d).^(1:K)) + d*((g1-d).^(1:K)+(g1+d).^(1:K)).*(1:K))/d^3;
%         v2 = 2.^(1-(1:K)).*((g1-d).^(1:K) - (g1+d).^(1:K) + d*((g1-d).^(0:K-1) + (g1+d).^(0:K-1)).*(1:K))/d^3;
%         G_g1 = toeplitz(sparse(1:K,1,v1(:),T,1),sparse(1,T));
%         G_g2 = toeplitz(sparse(1:K,1,v2(:),T,1),sparse(1,T));
%         grad = -((y-Gs)'*[G_g1*s,G_g2*s])';
%     end
% end