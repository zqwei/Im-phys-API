%
% this function runs our various oopsi filters, saves the results, and
% plots the inferred spike trains.  make sure that fast-oopsi and
% smc-oopsi repository are in your path if you intend to use them.
%
% to use the code, simply provide F, a vector of fluorescence observations,
% for each cell.  the fast-oopsi code can handle a matrix input,
% corresponding to a set of pixels, for each time bin. smc-oopsi expects a
% 1D fluorescence trace.
%
% see documentation for fast-oopsi and smc-oopsi to determine how to set
% variables
%
% input
%   F: fluorescence from a single neuron
%   fast_iter_max: max number of
%   P: Parameters of the model (optional)
%
% possible outputs
%   smc:    smc-oopsi estimate of {P[X_t|F]}_{t<T}, where X={n,C} or {n,C,h}, (smc.E), parameter estimates (smc.P), and structure of variables for algorithm (fast.V)
%
%
% Original by Josh Vogelstein
% https://github.com/jovo/smc-oopsi
% Revised by Ziqiang Wei
% weiz AT janelia DOT hhmi DOT org
%

%% main code
function [smc, fast, F, estF] = smc_oopsi(F, V, P)
if V.fast_iter_max < 1; V.fast_iter_max = 1; end;
if V.smc_iter_max < 1;  V.smc_iter_max = 1; end;
% preprocess - remove the lowest 10 frequencies
if V.preprocess==1
    V.T     = length(F);
    f       = detrend(F);
    nfft    = 2^nextpow2(V.T);
    y       = fft(f,nfft);
    bw      = 3;
    y(1:bw) = 0; y(end-bw+1:end)=0;
    iy      = ifft(y,nfft);
    F       = z1(real(iy(1:V.T)));
end
siz=size(F); if siz(1)>1, F=F'; end
% infer spikes using fast-oopsi (initializing pramaters for smc)
fprintf('\ninitializing pramaters for smc using fast-oopsi\n')
[fast.n, fast.P, fast.V, fast.C]= fast_oopsi(F,V);
% infer spikes using smc-oopsi (initializing pramaters for smc)
fprintf('starting smc-oopsi\n')
% bign1    = find(fast.n>0.1);
% bign1    = bign1(bign1>2);
% bign0    = bign1-1;
% df       = max((F(bign1)-F(bign0))./(F(bign0)));         % dff
P.C_0    = 0;                                            % baseline [Ca++] (\mu M)
P.C_init = P.C_0;                                        % initial [Ca++] (\mu M)
% S0       = Hill_v1(P,P.C_0);
% arg      = S0 + df*(S0 + 1/13);
% P.A      = ((arg*P.k_d)./(1-arg)).^(1/P.n)-P.C_0;
P.tau_c  = fast.V.dt/(1-fast.P.gam);                     % calcium decay time constant (sec)
nnorm    = V.n_max*fast.n/max(fast.n);                   % normalize inferred spike train
C        = filter(1,[1 -fast.P.gam],P.A*nnorm)'+P.C_0;   % calcium concentration
C1       = [Hill_v1(P,C); ones(1,V.T)];                  % for brevity
ab       = C1'\F';                                       % estimate scalse and offset
P.alpha  = ab(1);                                        % fluorescence scale: P.alpha = mean(F)
P.beta   = ab(2);                                        % fluorescence offset: P.beta  = min(F)
P.zeta   = (mad(F-ab'*C1,1)*1.4785)^2;
P.gamma  = P.zeta/5;                                     % signal dependent noise
% ZW: code from fast-oopsi git
k        = sum(nnorm);
P.k      = log(-log(1-k/V.T)/V.dt);             % baseline firing rate parameter
% ZW: what is spikegen.EFGinv???
% P.k      = V.spikegen.EFGinv(0.01, P, V);
[smc.M, smc.S, smc.P, smc.V] = smc_oopsi_v1(F,V,P);
S        = smc.S;
P        = smc.P;
smc.n_est= smc.M.nbar;
% ZW:infer the estimated fluroscence data
estF     = P.alpha*Hill_v1(P,sum(S.w_b.*S.C,1))+P.beta;

end

%% fast - oopsi
function [n_best, P_best, V_best, C_best]=fast_oopsi(F,V)
siz       = size(F);  if siz(2)==1, F=F'; end
% normalize F if it is only a trace
if V.Npixels==1
    F=detrend(F);
    F=F-min(F)+eps;
end
P.sig   = mean(mad(F',1)*1.4826);
P.gam   = (1-V.dt/1)*ones(V.Ncells,1);
P.lam   = 10*ones(V.Ncells,1);
P.a     = median(F,2);
if V.Npixels==1, P.b = quantile(F,0.05); else P.b=median(F,2); end
% for brevity and expediency
Z   = zeros(V.Ncells*V.T,1);                    % zero vector
M   = spdiags([repmat(-P.gam,V.T,1) repmat(Z,1,V.Ncells-1) (1+Z)], -V.Ncells:0,V.Ncells*V.T,V.Ncells*V.T);  % matrix transforming calcium into spikes, ie n=M*C
I   = speye(V.Ncells*V.T);                      % create out here cuz it must be reused
H1  = I;                                        % initialize memory for Hessian matrix
H2  = I;                                        % initialize memory for other part of Hessian matrix
d0  = 1:V.Ncells*V.T+1:(V.Ncells*V.T)^2;        % index of diagonal elements of TxT matrices
d1  = 1+V.Ncells:V.Ncells*V.T+1:(V.Ncells*V.T)*(V.Ncells*(V.T-1)); % index of off-diagonal elements of TxT matrices
posts = Z(1:V.fast_iter_max);                   % initialize likelihood
if numel(P.lam)==V.Ncells
    lam = V.dt*repmat(P.lam,V.T,1);             % for lik
elseif numel(P.lam)==V.Ncells*V.T
    lam = V.dt*P.lam;
else
    error('lam must either be length V.T or 1');
end

% if not iterating to estimate parameters, only this is necessary
[n, C, posts(1), V, M, d1, H1, d0, H2] = est_MAP(F,P, V, lam, M, d1, H1, d0, H2);
n_best = n;
P_best = P;
C_best = C;
V.fast_iter_tot = 1;
V.post = posts(1);
post_max = posts(1);

if V.fast_iter_max>1
    options = optimset('Display','off');        %#ok<NASGU> % don't show warnings for parameter estimation
    i       = 1;                                % iteration #
    % i_best  = i;                              %#ok<NASGU> % iteration with highest likelihood
    conv    = 0;                                % whether algorithm has converged yet
else
    conv    = 1;
end

V_best  = V;
while conv == 0
    i               = i+1;                      % update iteratation number
    V.fast_iter_tot = i;                        % record of total # of iterations    
    [P, V, lam, Z] = est_params(n, C, F, P, V, Z);% update parameters based on previous iteration
    [n, C, posts(i), V, M, d1, H1, d0, H2] = ...
        est_MAP(F,P, V, lam, M, d1, H1, d0, H2);% update inferred spike train based on new parameters
    
    if posts(i)>post_max || V.fast_ignore_post==1% if this is the best one, keep n and P
        n_best  = n;                            % keep n
        P_best  = P;                            % keep P
        C_best  = C;                            % keep C
        V_best  = V;
        % i_best  = i;                          %#ok<NASGU> % keep track of which was best
        post_max= posts(i);                     % keep max posterior
    end
    
    % if lik doesn't change much (relatively), or returns to some previous state, stop iterating
    if  i>=V.fast_iter_max || (abs((posts(i)-posts(i-1))/posts(i))<1e-3 || any(posts(1:i-1)-posts(i))<1e-5)
        % abs((posts(i)-posts(i-1))/posts(i))<1e-5 || posts(i-1)-posts(i)>1e5;
        disp('convergence criteria met')
        V.post  = posts(1:i);
        conv    = 1;
    end
end

V_best      = orderfields(V_best);                   % order fields alphabetically to they are easier to read
P_best      = orderfields(P_best);
% n_best      = n_best./repmat(max(n_best),V.T,1);
end

function [n, C, post, V, M, d1, H1, d0, H2] = est_MAP(F,P, V, lam, M, d1, H1, d0, H2)
% initialize n and C
z = 1;                                  % weight on barrier function
llam = reshape(1./lam',1,V.Ncells*V.T)';
n = 0.01+0*llam;                    % initialize spike train
C = 0*n;                                % initialize calcium
for j=1:V.Ncells
    C(j:V.Ncells:end) = filter(1,[1, -P.gam(j)],n(j:V.Ncells:end)); %(1-P.gam(j))*P.b(j);
end

% precompute parameters required for evaluating and maximizing likelihood
b       = repmat(P.b,1,V.T);           % for lik
M(d1)   = -repmat(P.gam,V.T-1,1);      % matrix transforming calcium into spikes, ie n=M*C
ba      = P.a'*b; ba=ba(:);            % for grad
aa      = repmat(diag(P.a'*P.a),V.T,1);% for grad
aF      = P.a'*F; aF=aF(:);            % for grad
e       = 1/(2*P.sig^2);               % scale of variance
H1(d0)  = -2*e*aa;                     % for Hess
grad_lnprior  = M'*llam;               % for grad


% find C = argmin_{C_z} lik + prior + barrier_z
while z>1e-13                           % this is an arbitrary threshold
    S = C;
    D = F-P.a*(reshape(S,V.Ncells,V.T))-b; % difference vector to be used in likelihood computation
    lik = e*D(:)'*D(:);                    % lik
    post = lik + llam'*n - z*sum(log(n));
    s    = 1;                                  % step size
    d    = 1;                                  % direction
    while norm(d)>5e-2 && s > 1e-3             % converge for this z (again, these thresholds are arbitrary)
        glik    = -2*e*(aF-aa.*C-ba);       % gradient
        g       = glik + grad_lnprior - z*M'*(n.^-1);
        H2(d0)  = n.^-2;                        % log barrier part of the Hessian
        H       = H1 - z*(M'*H2*M);             % Hessian
        d   = H\g;                              % direction to step using newton-raphson
        hit = -n./(M*d);                        % step within constraint boundaries
        hit=hit(hit>0);
        if any(hit<1)
            s = min(1,0.99*min(hit));
        else
            s = 1;
        end
        post1 = post+1;
        while post1>=post+1e-7                  % make sure newton step doesn't increase objective
            C1  = C+s*d;
            n   = M*C1;
            S1 = C1;
            D = F-P.a*(reshape(S1,V.Ncells,V.T))-b; % difference vector to be used in likelihood computation
            lik1 = e*D(:)'*D(:);             % lik
            post1 = lik1 + llam'*n - z*sum(log(n));
            s   = s/5;                  % if step increases objective function, decrease step size
            if s<1e-20; disp('reducing s further did not increase likelihood'), break; end      % if decreasing step size just doesn't do it
        end
        C    = C1;                      % update C
        post = post1;                   % update post
    end
    z=z/10;                             % reduce z (sequence of z reductions is arbitrary)
end

% reshape things in the case of multiple neurons within the ROI
n=reshape(n,V.Ncells,V.T)';
C=reshape(C,V.Ncells,V.T)';
end

function [P, V, lam, Z] = est_params(n, C, F, P, V, Z)
% generate regressor for spatial filter
% ZW: no nonlinear parameter or poisson parameter is handled here!
if V.fast_thr==1
    CC=0*C;
    for j=1:V.Ncells
        nsort   = sort(n(:,j));
        nthr    = nsort(round(0.98*V.T));
        nn      = Z(1:V.T);
        nn(n(:,j)<=nthr)=0;
        nn(n(:,j)>nthr)=1;
        CC(:,j) = filter(1,[1 -P.gam(j)],nn) + (1-P.gam(j))*P.b(j);
    end
else
    CC      = C;
end
A = [CC -1+Z(1:V.T)];
X = A\F';
% estimate - alpha
P.a = X(1:V.Ncells,:)';
% estimate - beta
P.b = X(end,:)';
b   = repmat(P.b,1,V.T);
% estimate - sigma
D   = F-P.a*(reshape(C,V.Ncells,V.T)) - b;
mse = D(:)'*D(:);
% estimate - sigma
P.sig = sqrt(mse)/V.T;
% estimate - lambda
nnorm   = n./repmat(max(n),V.T,1);
if numel(P.lam)==V.Ncells
    P.lam   = sum(nnorm)'/(V.T*V.dt);
    lam     = repmat(P.lam,V.T,1)*V.dt;
else
    P.lam   = nnorm/(V.T*V.dt);
    lam     = P.lam*V.dt;
end
end


%% smc - oopsi
function [M_best, S_best, P_best, V_best] = smc_oopsi_v1(F,V,P)
% initialize model Parameters
P.sigma_c       = 0.1;                                  % standard deviation of noise (\mu M)
if V.Nspikehist ==1                                     % if there are spike history terms
    P.omega     = -1;                                   % weight
    P.tau_h     = 0.02;                                 % time constant
    P.sigma_h   = 0;                                    % stan dev of noise
    P.g         = V.dt/P.tau_h;                         % for brevity
    P.sig2_h    = P.sigma_h^2*V.dt;                     % for brevity
end
P.a             = V.dt/P.tau_c;                         % for brevity
P.sig2_c        = P.sigma_c^2*V.dt;                     % for brevity
% initialize other stuff
P.lik       = -inf;                                     % we are trying to maximize the likelihood here
F           = max(F,eps);                               % in case there are any zeros in the F time series
S           = smc_oopsi_forward(F,V,P);                 % forward step
M           = smc_oopsi_backward(S,V,P);                % backward step
if V.smc_iter_max>1, P.conv=false; else P.conv=true; end
M_best      = M;
S_best      = S;
V_best      = V;
P_best      = P;


while P.conv==false
    post_max = max(P.lik);                        % store most recent parameter structure
    P        = update_params(V,S,M,P,F);          % update parameters + m step
    posts    = P.lik;
    i        = length(posts);
    S        = smc_oopsi_forward(F,V,P);          % forward step
    M        = smc_oopsi_backward(S,V,P);         % backward step
    
    if posts(i) > post_max                        % if this is the best one, keep n and P
        M_best  = M;                              % keep M
        S_best  = S;                              % keep S
        V_best  = V;                              % keep V
        P_best  = P;                              % keep P
    end
    
    if i>=V.smc_iter_max  || (abs((posts(i)-posts(i-1))/posts(i))<1e-3 || any(posts(1:i-1)-posts(i))<1e-5)        
        disp('convergence criteria met')
        P.conv=true;
    end
end
V.smc_iter_tot  = length(P.lik);
V_best          = orderfields(V_best);
P_best          = orderfields(P_best);
end

%% smc - forward
function S = smc_oopsi_forward(F,V,P)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that V.Nspikehist=0 (M is the # of spike history terms).
%
% The backward sampler has standard variance, as approximated typically.
% Each particle has the SAME backwards sampler, initialized at E[h_t]
% this function only does spike history stuff if necessary
%
% the model is F_t = f(C) = alpha C^n/(C^n + k_d) + beta + e_t,
% where e_t ~ N[0, gamma*f(C)+zeta]
%
% Inputs---
% F: Fluorescence
% V: Variables for algorithm to run
% P: initial Parameter estimates
%
% Outputs---
% S: simulation states
% allocate memory and initialize stuff
P.kx        = P.k'*V.x;
% extize particle info
S.p         = zeros(V.Nparticles,V.T);                  % extize rate
S.n         = false(V.Nparticles,V.T);                  % extize spike counts
S.C         = P.C_init*ones(V.Nparticles,V.T);          % extize calcium
S.w_f       = 1/V.Nparticles*ones(V.Nparticles,V.T);    % extize forward weights
S.w_b       = 1/V.Nparticles*ones(V.Nparticles,V.T);    % extize forward weights
S.Neff      = 1/V.Nparticles*ones(1,V.T_o);             % extize N_{eff}
% preprocess stuff for stratified resampling
ints        = linspace(0,1,V.Nparticles+1);             % generate intervals
diffs       = ints(2)-ints(1);                          % generate interval size
A.U_resamp  = repmat(ints(1:end-1),V.T_o,1)+diffs*rand(V.T_o,V.Nparticles); % resampling matrix
% extize misc stuff
A.U_sampl   = rand(V.Nparticles,V.T);                   % random samples
A.epsilon_c = sqrt(P.sig2_c)*randn(V.Nparticles,V.T);   % generate noise on C
% extize misc stuff
A.U_sampl   = rand(V.Nparticles,V.T);                   % random samples
A.epsilon_c = sqrt(P.sig2_c)*randn(V.Nparticles,V.T);   % generate noise on C

if V.Nspikehist>0                                   % if spike histories
    S.h = zeros(V.Nparticles,V.T,V.Nspikehist);     % extize spike history terms
    A.epsilon_h = zeros(V.Nparticles, V.T, V.Nspikehist); % generate noise on h
    for m=1:V.Nspikehist                            % add noise to each h
        A.epsilon_h(:,:,m) = sqrt(P.sig2_h(m))*randn(V.Nparticles,V.T);
    end
else                                                % if not, comput P[n_t] for all t
    S.p = repmat(1-exp(-exp(P.kx)*V.dt)',1,V.Nparticles)';
end

% extize stuff needed for conditional sampling
A.n_sampl   = rand(V.Nparticles,V.T);       % generate random number to use for sampling n_t
A.C_sampl   = rand(V.Nparticles,V.T);       % generate random number to use for sampling C_t
A.oney      = ones(V.Nparticles,1);         % vector of ones to call for various purposes to speed things up
A.zeroy     = zeros(V.Nparticles,1);        % vector of zeros

% extize stuff needed for REAL backwards sampling
O.p_o       = zeros(2^(V.freq-1),V.freq);   % extize backwards prob
O.mu_o      = zeros(2^(V.freq-1),V.freq);   % extize backwards mean
O.sig2_o    = zeros(1,V.freq);              % extize backwards variance
% extize stuff needed for APPROX backwards sampling
O.p         = zeros(V.freq,V.freq);         % extize backwards prob
O.mu        = zeros(V.freq,V.freq);         % extize backwards mean
O.sig2      = zeros(V.freq,V.freq);         % extize backwards var
% initialize backwards distributions
s           = V.freq;                       % initialize time of next observation
O.p_o(1,s)  = 1;                            % initialize P[F_s | C_s]
[O.mu_o(1,s), O.sig2_o(s)]  = init_lik(P,F(s));
O.p(1,s)    = 1;                            % initialize P[F_s | C_s]
O.mu(1,s)   = O.mu_o(1,s);                  % initialize mean of P[O_s | C_s]
O.sig2(1,s) = O.sig2_o(s);                  % initialize var of P[O_s | C_s]
if V.freq>1                                 % if intermitent sampling
    for tt=s:-1:2                           % generate spike binary matrix
        A.spikemat(:,tt-1) = repmat([zeros(1,2^(s-tt)) ones(1,2^(s-tt))],1,2^(tt-2))';
    end
    nspikes = sum(A.spikemat')';            %#ok<UDIM> % count number of spikes at each time step

    for n=0:V.freq-1
        A.ninds{n+1}= find(nspikes==n);     % get index for each number of spikes
        A.lenn(n+1) = length(A.ninds{n+1}); % find how many spikes
    end
end

O = update_moments(V,F,P,S,O,A,s);          % recurse back to get P[O_s | C_s] before the first observation

% do the particle filter
for t=V.freq+1:V.T-V.freq
    if V.condsamp==0 || (F(s+V.freq)-P.beta)/P.alpha>0.98 ||(F(s+V.freq)-P.beta)/P.alpha<0
        S = prior_sampler(V,F,P,S,A,t);     % use prior sampler when gaussian approximation is not good
    else                                    % otherwise use conditional sampler
        S = cond_sampler(V,F,P,S,O,A,t,s);
    end

    S.C(:,t) = S.next_C;                      % done to speed things up for older
    S.n(:,t) = S.next_n;                      % matlab, having issues with from-function
    if(isfield(S,'next_w_f'))               % update calls to large structures
        S.w_f(:,t)=S.next_w_f;
    else
        S.w_f(:,t)=S.w_f(:,t-1);
    end
    if V.Nspikehist>0                                % update S.h & S.p
        for m=1:V.Nspikehist, S.h(:,t,m)=S.h_new(:,1,m); end
        S.p(:,t)=S.p_new;
    end
    % at observations
    if mod(t,V.freq)==0
        % STRAT_RESAMPLE -- THIS WAS CAUSING TROUBLE IN OLDER MATLAB WITH
        % MANAGING LARGE ARRAYS IN S-STRUCTURE WITHIN THE FUNCTION CALL
        % S = strat_resample(V,S,t,A.U_resamp); % stratified resample
        Nresamp=t/V.freq;                                       % increase sample counter
        S.Neff(Nresamp)  = 1/sum(S.w_f(:,t).^2);                % store N_{eff}
        % if weights are degenerate or we are doing prior sampling then resample
        if S.Neff(Nresamp) < V.Nparticles/2 || V.condsamp==0
            [~, ind]    = histc(A.U_resamp(Nresamp,:),[0  cumsum(S.w_f(:,t))']);
            [~,ri]      = sort(rand(V.Nparticles,1));                    % these 3 lines stratified resample
            ind         = ind(ri);

            S.p(:,t-V.freq+1:t)   = S.p(ind,t-V.freq+1:t);      % resample probabilities (necessary?)
            S.n(:,t-V.freq+1:t)   = S.n(ind,t-V.freq+1:t);      % resample calcium
            S.C(:,t-V.freq+1:t)   = S.C(ind,t-V.freq+1:t);      % resample calcium
            S.w_f(:,t-V.freq+1:t) = 1/V.Nparticles*ones(V.Nparticles,V.freq); % reset weights
            if V.Nspikehist>0                                   % if spike history terms
                S.h(:,t-V.freq+1:t,:) = S.h(ind,t-V.freq+1:t,:);% resample all h's
            end
        end %function
        O = update_moments(V,F,P,S,O,A,t);                      % estimate P[O_s | C_tt] for all t'<tt<s as a gaussian
        s = t;                                                  % store time of last observation
    end
end %for time loop

end %conditional sampler

function [mu1, sig1] = init_lik(P,F)
% initialize likelihood
% get the mean (mu1) and variance (sig1) for P[C_t | F_t]
% compute mean
finv     = ((P.k_d*(F-P.beta))./(P.alpha-F+P.beta)).^(1/P.n); %initialize search with f^{-1}(o)
mu1      = finv;
if mu1>0 && imag(mu1)==0
    sig1 = -1/(-(-P.alpha*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)+P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)+2*(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2*(-P.alpha*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)+P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)-(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)*(-P.alpha*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)+P.alpha*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)+3*P.alpha*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2-2*P.alpha*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2-P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)-(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^3*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2+1/2*(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2*(P.gamma*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)-P.gamma*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)-3*P.gamma*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2+2*P.gamma*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2+P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)-1/2*(P.gamma*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)-P.gamma*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)-3*P.gamma*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2+2*P.gamma*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2+P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)+1/2*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2);
else
    mu1  = 0;
    sig1 = 0;
end
end %init_lik

function O = update_moments(V,F,P,S,O,A,t)
% update moments
%%%% maybe make a better proposal for epi
s           = V.freq;                   % find next observation time
[mu1, sig1] = init_lik(P,F(t+s));
O.mu_o(1,s) = mu1;                      % initialize mean of P[O_s | C_s]
O.sig2_o(s) = sig1;                     % initialize var of P[O_s | C_s]
O.p(1,s)    = 1;                        % initialize P[F_s | C_s]
O.mu(1,s)   = mu1;                      % initialize mean of P[O_s | C_s]
O.sig2(1,s) = sig1;                     % initialize var of P[O_s | C_s]
if V.Nspikehist>0
    hhat        = zeros(V.freq,V.Nspikehist);                   % extize hhat
    phat        = zeros(1,V.freq+1);                            % extize phat
    hs          = S.h(:,t,:);                                   % this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:V.Nspikehist)= hs(:,1,1:V.Nspikehist);                % this too
    hhat(1,:)   = sum(repmat(S.w_f(:,t),1,V.Nspikehist).*h,1);  % initialize hhat
    phat(1)     = sum(S.w_f(:,t).*S.p(:,t),1);                  % initialize phat
end
if V.Nspikehist>0
    for tt=1:s
        % update hhat
        for m=1:V.Nspikehist                                    % for each spike history term
            hhat(tt+1,m)=(1-P.g(m))*hhat(tt,m)+phat(tt);
        end
        y_t         = P.kx(tt+t)+P.omega'*hhat(tt+1,:)';        % input to neuron
        phat(tt+1)  = 1-exp(-exp(y_t)*V.dt);                    % update phat
    end
else
    phat  = 1-exp(-exp(P.kx(t+1:t+s)')*V.dt);                   % update phat
end
for tt=s:-1:2
    O.p_o(1:2^(s-tt+1),tt-1)    = repmat(O.p_o(1:2^(s-tt),tt),2,1).*[(1-phat(tt))*ones(1,2^(s-tt)) phat(tt)*ones(1,2^(s-tt))]';
    O.mu_o(1:2^(s-tt+1),tt-1)   = (1-P.a)^(-1)*(repmat(O.mu_o(1:2^(s-tt),tt),2,1)-P.A*A.spikemat(1:2^(s-tt+1),tt-1)-P.a*P.C_0);     %mean of P[O_s | C_k]
    O.sig2_o(tt-1)              = (1-P.a)^(-2)*(P.sig2_c+O.sig2_o(tt)); % var of P[O_s | C_k]
    for n=0:s-tt+1
        nind=A.ninds{n+1};
        O.p(n+1,tt-1)   = sum(O.p_o(nind,tt-1));
        ps          = (O.p_o(nind,tt-1)/O.p(n+1,tt-1))';
        O.mu(n+1,tt-1)  = ps*O.mu_o(nind,tt-1);
        O.sig2(n+1,tt-1)= O.sig2_o(tt-1) + ps*(O.mu_o(nind,tt-1)-repmat(O.mu(n+1,tt-1)',A.lenn(n+1),1)).^2;
    end
end
if s==2
    O.p     = O.p_o;
    O.mu    = O.mu_o;
    O.sig2  = repmat(O.sig2_o,2,1);
    O.sig2(2,2) = 0;
end
while any(isnan(O.mu(:)))              % in case ps=0/0, which yields a NaN, approximate mu and sig
    O.mu(1,:)   = O.mu_o(1,:);
    O.sig2(1,:) = O.sig2_o(1,:);
    ind         = find(isnan(O.mu));
    O.mu(ind)   = O.mu(ind-1)-P.A;
    O.sig2(ind) = O.sig2(ind-1);
end
O.p=O.p+eps; % such that there are no actual zeros
end %function UpdateMoments

function S = prior_sampler(V,F,P,S,A,t)
% particle filtering using the prior sampler
if V.Nspikehist>0                                % update noise on h
    S.h_new=zeros(size(S.n,1),1,V.Nspikehist);
    for m=1:V.Nspikehist
        S.h_new(:,1,m)=(1-P.g(m))*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end
    % update rate and sample spikes
    hs          = S.h_new;              % this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:V.Nspikehist)  = hs(:,1,1:V.Nspikehist);        % this too
    y_t         = P.kx(t)+P.omega'*h';  % input to neuron
    S.p_new     = 1-exp(-exp(y_t)*V.dt);% update rate for those particles with y_t<0
    S.p_new     = S.p_new(:);
else
    S.p_new     = S.p(:,t);
end
if ~V.use_true_n
    S.next_n    = A.U_sampl(:,t)<S.p_new;   % sample n
else
    S.next_n    = V.true_n;
end
S.next_C        = (1-P.a)*S.C(:,t-1)+P.A*S.next_n+P.a*P.C_0+A.epsilon_c(:,t);% sample C
% get weights at every observation          %THIS NEEDS FIX FOR EPI DATA
if mod(t,V.freq)==0
    S_mu        = Hill_v1(P,S.next_C);
    F_mu        = P.alpha*S_mu+P.beta;      % compute E[F_t]
    F_var       = P.gamma*S_mu+P.zeta;      % compute V[F_t]
    %%%% this must also change for epi
    ln_w        = -0.5*(F(t)-F_mu).^2./F_var - log(F_var)/2;% compute log of weights
    ln_w        = ln_w-max(ln_w);           % subtract the max to avoid rounding errors
    w           = exp(ln_w);                % exponentiate to get actual weights
    %     error('forgot to include the previous weight in this code!!!!')
    %     break
    S.next_w_f  = w/sum(w);                 % normalize to define a legitimate distribution
end
end

function S = cond_sampler(V,F,P,S,O,A,t,s)
% particle filtering using the CONDITIONAL sampler
% if spike histories, sample h and update p
if V.Nspikehist>0                                    % update noise on h
    S.h_new=zeros(size(S.n,1),1,V.Nspikehist);
    for m=1:V.Nspikehist                             % for each spike history term
        S.h_new(:,1,m)=(1-P.g(m))*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end
    hs              = S.h_new;              % this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:V.Nspikehist)      = hs(:,1,1:V.Nspikehist);        % this too
    S.p_new         = 1-exp(-exp(P.kx(t)+P.omega'*h')*V.dt);% update p
    S.p_new         = S.p_new(:);
else
    S.p_new         = S.p(:,t);
end
% compute P[n_k | h_k]
ln_n    = [log(S.p_new) log(1-S.p_new)];    % compute [log(spike) log(no spike)]
% compute log G_n(n_k | O_s) for n_k=1 and n_k=0
k   = V.freq-(t-s)+1;
m0  = (1-P.a)*S.C(:,t-1)+P.a*P.C_0;         % mean of P[C_k | C_{t-1}, n_k=0]
m1  = (1-P.a)*S.C(:,t-1)+P.A+P.a*P.C_0;     % mean of P[C_k | C_{t-1}, n_k=1]
m2  = O.mu(1:k,t-s);                        % mean of P[O_s | C_k] for n_k=1 and n_k=0
v2  = O.sig2(1:k,t-s);                      % var of P[O_s | C_k] for n_k=1 and n_k=0
v   = repmat(P.sig2_c+v2',V.Nparticles,1);           % var of G_n(n_k | O_s) for n_k=1 and n_k=0
ln_G0= -0.5*log(2*pi.*v)-.5*(repmat(m0,1,k)-repmat(m2',V.Nparticles,1)).^2./v;   % log G_n(n_k | O_s) for n_k=1 and n_k=0
ln_G1= -0.5*log(2*pi.*v)-.5*(repmat(m1,1,k)-repmat(m2',V.Nparticles,1)).^2./v;   % log G_n(n_k | O_s) for n_k=1 and n_k=0
mx  = max(max(ln_G0,[],2),max(ln_G1,[],2))';% get max of these
mx  = repmat(mx,k,1)';
G1  = exp(ln_G1-mx);                        % norm dist'n for n=1;
M1  = G1*O.p(1:k,t-s);                      % times prob of n=1
G0  = exp(ln_G0-mx);                        % norm dist'n for n=0;
M0  = G0*O.p(1:k,t-s);                      % times prob n=0
ln_G    = [log(M1) log(M0)];                % ok, now we actually have the gaussians
% compute q(n_k | h_k, O_s)
ln_q_n  = ln_n + ln_G;                      % log of sampling dist
mx      = max(ln_q_n,[],2);                 % find max of each column
mx2     = repmat(mx,1,2);                   % matricize
q_n     = exp(ln_q_n-mx2);                  % subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
q_n     = q_n./repmat(sum(q_n,2),1,2);      % normalize to make this a true sampling distribution (ie, sums to 1)
% sample n
S.next_n= A.n_sampl(:,t)<q_n(:,1);          % sample n
sp      = S.next_n==1;                      % store index of which samples spiked
nosp    = S.next_n==0;                      % and which did not
% sample C
if mod(t,V.freq)==0                         % if not intermittent
    v       = repmat(O.sig2(1,t-s),V.Nparticles,1);  % get var
    m       = repmat(O.mu(1,t-s),V.Nparticles,1);    % get mean
else                                        % if intermittent, sample from mixture
    % first sample component
    if(isempty(find(sp,1))), sp_i=[];       % handler for empty spike trains
    else [~,sp_i]   = histc(A.C_sampl(sp,t),[0  cumsum(O.p(1:k-1,t-s))'/sum(O.p(1:k-1,t-s))]); end
    if(isempty(find(nosp,1))), nosp_i=[];   % handle for saturated spike trains
    else [~,nosp_i] = histc(A.C_sampl(nosp,t),[0  cumsum(O.p(1:k,t-s))'/sum(O.p(1:k,t-s))]); end
    v       = O.sig2(1:k,t-s);              % get var of each component
    v(sp)   = v(sp_i);                      % if particle spiked, then use variance of spiking
    v(nosp) = v(nosp_i);                    % o.w., use non-spiking variance
    m       = O.mu(1:k,t-s);                % get mean of each component
    m(sp)   = m(sp_i);                      % if particle spiked, use spiking mean
    m(nosp) = m(nosp_i);                    % o.w., use non-spiking mean
end
v_c         = (1./v+1/P.sig2_c).^(-1);      %variance of dist'n for sampling C
m_c         = v_c.*(m./v+((1-P.a)*S.C(:,t-1)+P.A*S.next_n+P.a*P.C_0)/P.sig2_c);%mean of dist'n for sampling C
S.next_C    = normrnd(m_c,sqrt(v_c));       % sample C
% update weights
if mod(t,V.freq)==0                         % at observations compute P(O|H)
    S_mu        = Hill_v1(P,S.next_C);
    if V.scan==0                            % when doing epi, also add previous time steps
        for tt=s+1:t-1, S_mu=S_mu+Hill_v1(P,S.C(s+tt)); end
    end
    F_mu        = P.alpha*S_mu+P.beta;      % compute E[F_t]
    F_var       = P.gamma*S_mu+P.zeta;      % compute V[F_t]
    %%%% log_PO_H must change for epi
    log_PO_H    = -0.5*(F(t)-F_mu).^2./F_var - log(F_var)/2; % compute log of weights
else                                        % when no observations are present, P[O|H^{(i)}] are all equal
    log_PO_H    = (1/V.Nparticles)*A.oney;
end
log_n           = A.oney;                   % extize log sampling spikes
log_n(sp)       = log(S.p_new(sp));         % compute log P(spike)
log_n(nosp)     = log(1-S.p_new(nosp));     % compute log P(no spike)
log_C_Cn        = -0.5*(S.next_C-((1-P.a)*S.C(:,t-1)+P.A*S.next_n+P.a*P.C_0)).^2/P.sig2_c;%log P[C_k | C_{t-1}, n_k]

log_q_n         = A.oney;                   % initialize log q_n
log_q_n(sp)     = log(q_n(sp,1));           % compute what was the log prob of sampling a spike
log_q_n(nosp)   = log(1-q_n(nosp,1));       % or sampling no spike
log_q_C         = -0.5*(S.next_C-m_c).^2./v_c;% log prob of sampling the C_k that was sampled

log_quotient    = log_PO_H + log_n + log_C_Cn - log_q_n - log_q_C;

sum_logs        = log_quotient+log(S.w_f(:,t-1));   % update log(weights)
w               = exp(sum_logs-max(sum_logs));      % exponentiate log(weights)
S.next_w_f      = w./sum(w);                        % normalize such that they sum to unity
if any(isnan(w))
    warning('smc:weights','some weights are NaN')
%     keyboard,
end
end %condtional sampler

%% smc - backward
function M = smc_oopsi_backward(S,V,P)
% this function iterates backward one step computing P[H_t | H_{t+1},O_{0:T}]
% Input---
% Sim:  simulation metadata
% S:    particle positions and weights
% P:    parameters
% Z:    a bunch of stuff initialized for speed
% t:    current time step
%
% Output is a single structure Z with the following fields
% n1:   vector of spikes or no spike for each particle at time t
% C0:   calcium positions at t-1
% C1:   calcium positions at t (technically, this need not be output)
% C1mat:matrix from C1
% C0mat:matrix from C0
% w_b:  backwards weights
Z.oney  = ones(V.Nparticles,1);                 % initialize stuff for speed
Z.zeroy = zeros(V.Nparticles);
Z.C0    = S.C(:,V.T);
Z.C0mat = Z.C0(:,Z.oney)';
if V.est_c==false                               % if not maximizing the calcium parameters, then the backward step is simple
    if V.use_true_n                             % when spike train is provided, backwards is not necessary
        S.w_b=S.w_f;
    else
        for t=V.T-V.freq-1:-1:V.freq+1          % actually recurse backwards for each time step
            Z = step_backward(V,S,P,Z,t);
            S.w_b(:,t-1) = Z.w_b;               % update forward-backward weights
        end
    end
else                                            % if maximizing calcium parameters,
    % need to compute some sufficient statistics
    M.Q = zeros(3);                             % the quadratic term for the calcium par
    M.L = zeros(3,1);                           % the linear term for the calcium par
    M.J = 0;                                    % remaining terms for calcium par
    M.K = 0;
    for t=V.T-V.freq-1:-1:V.freq+1
        if V.use_true_n                         % force true spikes hack
            Z.C0    = S.C(t-1);
            Z.C0mat = Z.C0;
            Z.C1    = S.C(t);
            Z.C1mat = Z.C1;
            Z.PHH   = 1;
            Z.w_b   = 1;
            Z.n1    = S.n(t);
        else
            Z = step_backward(V,S,P,Z,t);
        end
        S.w_b(:,t-1) = Z.w_b;
        % below is code to quickly get sufficient statistics
        C0dt    = Z.C0*V.dt;
        bmat    = Z.C1mat-Z.C0mat';
        bPHH    = Z.PHH.*bmat;
        M.Q(1,1)= M.Q(1,1) + sum(Z.PHH*(C0dt.^2));  % Q-term in QP
        M.Q(1,2)= M.Q(1,2) - Z.n1'*Z.PHH*C0dt;
        M.Q(1,3)= M.Q(1,3) + sum(sum(-Z.PHH.*Z.C0mat'*V.dt^2));
        M.Q(2,2)= M.Q(2,2) + sum(Z.PHH'*(Z.n1.^2));
        M.Q(2,3)= M.Q(2,3) + sum(sum(Z.PHH(:).*repmat(Z.n1,V.Nparticles,1))*V.dt);
        M.Q(3,3)= M.Q(3,3) + sum(Z.PHH(:))*V.dt^2;
        M.L(1)  = M.L(1) + sum(bPHH*C0dt);          % L-term in QP
        M.L(2)  = M.L(2) - sum(bPHH'*Z.n1);
        M.L(3)  = M.L(3) - V.dt*sum(bPHH(:));
        M.J     = M.J + sum(Z.PHH(:));              % J-term in QP /sum J^(i,j)_{t,t-1}/
        M.K     = M.K + sum(Z.PHH(:).*bmat(:).^2);  % K-term in QP /sum J^(i,j)_{t,t-1} (d^(i,j)_t)^2/
    end
    M.Q(2,1) = M.Q(1,2);                            % symmetrize Q
    M.Q(3,1) = M.Q(1,3);
    M.Q(3,2) = M.Q(2,3);
end
% copy particle swarm for later
M.w = S.w_b;
M.n = S.n;
M.C = S.C;
if isfield(S,'h'), M.h=S.h; end
M.nbar = sum(S.w_b.*S.n,1);
end

function Z = step_backward(V,S,P,Z,t)
% compute ln P[n_t^i | h_t^i]
Z.n1            = S.n(:,t);                         % for prettiness sake
ln_Pn           = 0*Z.oney;                         % for fastiness sake
ln_Pn(Z.n1==1)  = log(S.p(Z.n1==1,t));              % P[n=1] for those that spiked
ln_Pn(~Z.n1)    = log(1-S.p(~Z.n1,t));              % P[n=0] for those that did not
% compute ln P[C_t^i | C_{t-1}^j, n_t^i]
Z.C0        = S.C(:,t-1);                           % for prettiness sake
Z.C1        = S.C(:,t);
Z.C1mat     = Z.C1(:,Z.oney);                       % recall from previous time step
Z.C0mat     = Z.C0(:,Z.oney);                       % faster than repamt
mu          = (1-P.a)*S.C(:,t-1)+P.A*Z.n1+P.a*P.C_0;% mean
mumat       = mu(:,Z.oney)';                        % faster than repmat
ln_PC_Cn    = -0.5*(Z.C1mat - mumat).^2/P.sig2_c;   % P[C_t^i | C_{t-1}^j, n_t^i]
% compute ln P[h_t^i | h_{t-1}^j, n_{t-1}^i]
ln_Ph_hn    = Z.zeroy;                              % reset transition prob for h terms
for m=1:V.Nspikehist                                % for each h term
    h1      = S.h(:,t,m);                           % faster than repmat
    h1      = h1(:,Z.oney);
    h0      = P.g(m)*S.h(:,t-1,m)+S.n(:,t-1);
    h0      = h0(:,Z.oney)';
    ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/P.sig2_h(m);
end
% compute P[H_t^i | H_{t-1}^j]
sum_lns = ln_Pn(:,Z.oney)+ln_PC_Cn + ln_Ph_hn;      % in order to ensure this product doesn't have numerical errors
mx      = max(sum_lns,[],1);                        % find max in each of row
mx      = mx(Z.oney,:);                             % make a matrix of maxes
T0      = exp(sum_lns-mx);                          % exponentiate subtracting maxes (so that in each row, the max entry is exp(0)=1
Tn      = sum(T0,1);                                % then normalize
T       = T0.*repmat(1./Tn(:)', V.Nparticles, 1);   % such that each column sums to 1
% compute P[H_t^i, H_{t-1}^j | O]
PHHn    = (T*S.w_f(:,t-1))';                        % denominator
PHHn(PHHn==0) = eps;
PHHn2   = PHHn(Z.oney,:)';                          % faster than repmat
PHH     = T .* (S.w_b(:,t)*S.w_f(:,t-1)')./PHHn2;   % normalize such that sum(PHH)=1
sumPHH  = sum(PHH(:));
if sumPHH==0
    Z.PHH = ones(V.Nparticles)/(V.Nparticles);
else
    Z.PHH   =  PHH/sum(PHH(:));
end
Z.w_b   = sum(Z.PHH,1);                             % marginalize to get P[H_t^i | O]

if any(isnan(Z.w_b))
    return
end
end

%% other functions
function F = Hill_v1(P,C)
C(C<0)  = 0;
F       = C.^P.n./(C.^P.n+P.k_d);
end

function x = z1(y)

x = (y-min(y(:)))/(max(y(:))-min(y(:)))+eps;
end

function Enew = update_params(V,S,M,E,F)
    Enew       = E;    % initialize parameters
    lik        = [];   % initialize likelihood
    optionsQP  = optimset('Display','off');
    optionsGLM = optimset('Display','off');%, 'TolFun',1e-6, 'Algorithm', 'interior-point', 'Hessian', 'bfgs'); 
    % ZW: interior-point does not work well
    % ZW: using quasi-newton instead --> sqp
    % ZW: remve 'GradObj','off',

    % MLE for spiking parameters
    if V.est_n == true
        % MLE for spike rate parameters: baseline (b), linear filter (k), and spike history weights (omega)
        RateParams = E.k;                      % vector of parameters to estimate (changes depending on user input of which parameters to estimate)
        sp         = S.n==1;                   % find (particles,time step) pairs that spike
        nosp       = S.n==0;                   % don't spike
        x        = repmat(V.x,1,V.Nparticles); % generate matrix for gradinent
        zeroy      = zeros(V.Nparticles,V.T);  % make matrix of zeros for evaluating lik

        if V.est_h == true
            if V.Nspikehist>0                     % if spike history terms are present
                RateParams=[RateParams; E.omega]; % also estimate omega
                    for i=1:V.Nspikehist          % and modify stimulus matrix for gradient
                        x(V.StimDim+i,:)=reshape(S.h(:,:,i),1,V.Nparticles*V.T);
                    end
            end
            Z=ones(size(RateParams));
            [bko, lik_r]=fmincon(@(RateParams) f_bko(RateParams, V, S, sp, nosp, zeroy, x),...
                            RateParams,[],[],[],[],-10*Z,10*Z,[],optionsGLM); %fix for h-problem
            Enew.k      = bko(1:end-V.Nspikehist);                            % set new parameter estimes
            if V.Nspikehist>0, Enew.omega = bko(end-V.Nspikehist+1:end); end  % for omega too
        else
            if V.Nspikehist>0                   % if spike history terms are present
                for i=1:V.Nspikehist            % and modify stimulus matrix for gradient
                    x(V.StimDim+i,:)=reshape(S.h(:,:,i),1,V.Nparticles*V.T);
                end
            end

            [bk, lik_r]= fminunc(@(RateParams) f_bk(RateParams, V, S, sp, nosp, zeroy, x),...
                            RateParams, optionsGLM); % find MLE
            Enew.k     = bk(1:end);             % set new parameter estimes
        end
        Enew.lik_r   = -lik_r;
        lik = [lik Enew.lik_r];
    end
    
    % MLE for calcium parameters
    if V.est_c == true
        if V.est_t == 0
            [ve_x, fval] = quadprog(M.Q(2:3,2:3), M.L(2:3),[],[],[],[],[0 0],[inf inf],[E.A E.C_0/E.tau_c]+eps,optionsQP);
            Enew.tau_c  = E.tau_c;
            Enew.A      = ve_x(1);
            Enew.C_0    = ve_x(2)/E.tau_c;
        else
            [ve_x, fval] = quadprog(M.Q, M.L,[],[],[],[],[0 0 0],[inf inf inf],[1/E.tau_c E.A E.C_0/E.tau_c]+eps,optionsQP);
            Enew.tau_c  = 1/ve_x(1);
            Enew.A      = ve_x(2);
            Enew.C_0    = ve_x(3)/ve_x(1);
        end
        fval        = M.K/2 + fval;                         % variance
        Enew.sigma_c= sqrt(fval/(M.J*V.dt));                % factor in dt
        Enew.lik_c  = - fval/(Enew.sigma_c*sqrt(V.dt)) - M.J*log(Enew.sigma_c);
        lik = [lik Enew.lik_c];

        Enew.a         = V.dt/E.tau_c;                      % for brevity
        Enew.sig2_c    = E.sigma_c^2*V.dt;                  % for brevity
    end

    % MLE for spike history parameters
    for m=1:V.Nspikehist
        Enew.sigma_h(m)= sum(M.v{m})/V.T;
    end

    % MLE for observation parameters

    if V.est_F == true
        [Enew.lik_o, ab] = f1_ab(E, S, F,optionsQP);
        Enew.alpha       = ab(1);
        Enew.beta        = ab(2);
        Enew.zeta        = E.zeta*ab(3);
        Enew.gamma       = E.gamma*ab(3);
        lik              = [lik Enew.lik_o];
    end
    Enew.lik=[E.lik sum(lik)];
end

function [lik, dlik]= f_bko(RateParams, V, S, sp, nosp, zeroy, x)    % get lik and grad
xk      = RateParams(1:end-V.Nspikehist)'*V.x;      % filtered stimulus
hs      = zeroy;                                    % incorporate spike history terms
for l=1:V.Nspikehist
    hs  = hs+RateParams(end-V.Nspikehist+l)*S.h(:,:,l); 
end
s       = repmat(xk,V.Nparticles,1) + hs;
f_kdt   = exp(s)*V.dt;                              % shorthand
ef      = exp(f_kdt);                               % shorthand
lik     = -sum(S.w_b(sp).*log(1-1./ef(sp)))...      % liklihood
    +sum(S.w_b(nosp).*f_kdt(nosp));
if nargout > 1                                      % if gradobj=on
    dlik      = -x(:,sp)*(S.w_b(sp).*f_kdt(sp)./(ef(sp)-1))... %gradient of lik
        + x(:,nosp)*(S.w_b(nosp).*f_kdt(nosp));
end
end %function f_bko

function [lik, dlik]= f_bk(RateParams, V, S, sp, nosp, zeroy, x)              % get lik and grad
xk      = RateParams'*V.x;                          % filtered stimulus
hs      = zeroy;                                    % incorporate spike history terms
for l=1:V.Nspikehist, hs  = hs+E.omega*S.h(:,:,l); end
s       = repmat(xk,V.Nparticles,1) + hs;
f_kdt   = exp(s)*V.dt;                              % shorthand
ef      = exp(f_kdt);                               % shorthand
lik     = -sum(S.w_b(sp).*log(1-1./ef(sp)))...      % liklihood
    +sum(S.w_b(nosp).*f_kdt(nosp));
if nargout > 1                                      % if gradobj=on
    dlik      = -x(:,sp)*(S.w_b(sp).*f_kdt(sp)./( ef(sp)-1))... % gradient of lik
        + x(:,nosp)*(S.w_b(nosp).*f_kdt(nosp));
end
end %function f_bko

function [lik, x] = f1_ab(E, S, F,optionsQP)
%find MLE for {alpha, beta and gamma/zeta}
%THIS EXPLICITLY ASSUMES WEIGHTS w_b ARE SUM=1 NORMALIZED (!)
pfS         = Hill_v1(E,S.C);
pfV         = E.gamma*pfS+E.zeta;
% minimize quadratic form of E[(F - ab(1)*pfS - ab(2))^2/pfV]
% taken as weighted average over all particles (Fmean)
f1_abn      = sum(sum(S.w_b,2));       % normalization
f1_abH(1,1) = sum(sum(pfS.^2./pfV.*S.w_b,2));% QF
f1_abH(2,2) = sum(sum(S.w_b./pfV,2));
f1_abH(1,2) = sum(sum(pfS./pfV.*S.w_b,2));
f1_abH(2,1) = f1_abH(1,2);
f1_abf      = [0;0]; f1_abc=0;          % LF and offset
for i       = 1:size(pfS,1)             % over particles
    f1_abf(1)   = f1_abf(1) - sum(F(:)'.*pfS(i,:)./pfV(i,:).*S.w_b(i,:));
    f1_abf(2)   = f1_abf(2) - sum(F(:)'./pfV(i,:).*S.w_b(i,:));    
    f1_abc      = f1_abc + sum(F(:)'.^2./pfV(i,:).*S.w_b(i,:));
end
% solve QP given ab(1)>0, no bound on ab(2)
[x, lik]    = quadprog(f1_abH,f1_abf,[],[],[],[],[0 -inf],[inf inf],[],optionsQP);
lik         = (lik+f1_abc/2);             % variance
x(3)        = lik/f1_abn;                % estimate gamma_new/gamma
lik         = -lik-sum(sum(log(x(3)*pfV).*S.w_b,2))/2;
end %function f_ab
