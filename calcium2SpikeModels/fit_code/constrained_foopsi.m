function [c, b, c1, g, sn, sp] = constrained_foopsi(y,b,c1,g,sn,options)
% spike inference using a constrained foopsi approach:
%      min      sum(sp)
%    c,sp,b,c1
%      subject to: sp >= 0
%                   b >= 0
%                  G*c = sp
%                   c1 >= 0
%           ||y-b-c - c_in|| <= sn*sqrt(T)
%   Variables:
%   y:      raw fluorescence data (vector of length(T))
%   c:      denoised calcium concentration (Tx1 vector)
%   b:      baseline concentration (scalar)
%  c1:      initial concentration (scalar)
%   g:      discrete time constant(s) (scalar or 2x1 vector)
%  sn:      noise standard deviation (scalar)
%  sp:      spike vector (Tx1 vector)
%   USAGE:
%   [c,b,c1,g,sn,sp] = constrained_foopsi(y,b,c1,g,sn,OPTIONS)
%   The parameters b,cin,g,sn can be given or else are estimated from the data
%   OPTIONS: (stuct for specifying options)
%         p: order for AR model, used when g is not given (default 2)
%    method: methods for performing spike inference
%   available methods: 'dual' uses dual ascent
%                       'cvx' uses the cvx package available from cvxr.com (default)
%                      'lars' uses the least regression algorithm 
%                     'spgl1' uses the spgl1 package available from
%                     math.ucdavis.edu/~mpf/spgl1/  (usually fastest)
%   bas_nonneg:   flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y)
%   noise_range:  frequency range over which the noise power is estimated. Default [Fs/4,Fs/2]
%   noise_method: method to average the PSD in order to obtain a robust noise level estimate
%   lags:         number of extra autocovariance lags to be considered when estimating the time constants
%   resparse:     number of times that the solution is resparsened (default 0). Currently available only with methods 'cvx', 'spgl'
%   fudge_factor: scaling constant to reduce bias in the time constant estimation (default 1 - no scaling)
% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015 

defoptions.p = 2;
defoptions.method = 'cvx';
defoptions.bas_nonneg = 1;              % nonnegativity option for baseline estimation
defoptions.noise_range = [0.25,0.5];    % frequency range over which to estimate the noise
defoptions.noise_method = 'logmexp';    % method for which to estimate the noise level
defoptions.lags = 5;                    % number of extra lags when computing the AR coefficients
defoptions.resparse = 0;                % number of times to re-sparse solution
defoptions.fudge_factor = 1;            % fudge factor for time constants

if nargin < 6
    options = defoptions;
    if nargin < 5
        sn = [];
        if nargin < 4
            g = [];
            if nargin < 3
                c1 = [];
                if nargin < 2
                    b = [];
                end
            end
        end
    end
end
      
if ~isfield(options,'p');  options.p = defoptions.p;  end
if ~isfield(options,'method'); options.method = defoptions.method; end
if ~isfield(options,'bas_nonneg'); options.bas_nonneg = defoptions.bas_nonneg; end
if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'lags'); options.lags = defoptions.lags; end
if ~isfield(options,'resparse'); options.resparse = defoptions.resparse; end
if ~isfield(options,'fudge_factor'); options.fudge_factor = defoptions.fudge_factor; end

method = options.method;    
if isempty(b);
    bas_est = 1;
else
    bas_est = 0;
end
if isempty(c1)
    c1_est = 1;
else
    c1_est = 0;
end

y = y(:);
T = length(y);
y_full = y;
mis_data = isnan(y);
E = speye(T);
E(mis_data,:) = [];

if any(mis_data)
    y_full(mis_data) = interp1(find(~mis_data),y(~mis_data),find(mis_data));
end
    

if isempty(sn)
    sn = GetSn(y_full,options.noise_range,options.noise_method);
end
if isempty(g)
    g = estimate_time_constants(y_full,options.p,sn,options.lags);
    while max(abs(roots([1,-g(:)']))>1) && options.p < 5
        warning('No stable AR(%i) model found. Checking for AR(%i) model \n',options.p,options.p+1);
        options.p = options.p + 1;
        g = estimate_time_constants(y,options.p,sn,options.lags);
    end
    if options.p == 5
        g = 0;
    end
    %fprintf('Stable AR(%i) model found \n',options.p);
    % re-adjust time constant values
    rg = roots([1;-g(:)]);
    if ~isreal(rg); rg = real(rg) + .001*randn(size(rg)); end
    rg(rg>1) = 0.95 + 0.001*randn(size(rg(rg>1)));
    rg(rg<0) = 0.15 + 0.001*randn(size(rg(rg<0)));
    pg = poly(options.fudge_factor*rg);
    g = -pg(2:end);
end
if options.bas_nonneg  % lower bound for baseline
    b_lb = 0;
else
    b_lb = min(y);
end

if strcmpi(method,'dual'); method = 'dual';
elseif strcmpi(method,'cvx'); method = 'cvx';
elseif strcmpi(method,'lars'); method = 'lars';
elseif strcmpi(method,'spgl1'); method = 'spgl1';
else fprintf('Invalid choice of method. Using CVX \n'); method = 'cvx';
end

if strcmpi(method,'dual') && any(mis_data)
    warning('Dual method does not support missing data. Switching to CVX');
    method = 'cvx';
end

if options.resparse > 0 && (strcmpi(method,'dual') || strcmpi(method,'lars'))
    warning('Resparsening is not supported with chosen method. Switching to CVX');
    method = 'cvx';
end

pathCell = regexp(path, pathsep, 'split'); %#ok<NASGU>
g = g(:);
G = spdiags(ones(T,1)*[-g(end:-1:1)',1],-length(g):0,T,T);
gd = max(roots([1,-g']));  % decay time constant for initial concentration
gd_vec = gd.^((0:T-1)');

switch method
    case 'dual'
        v = G'*ones(T,1); %#ok<NASGU>
        thr = sn*sqrt(T-sum(mis_data));
        if bas_est; b = 0; end
        if c1_est; c1 = 0; end
        c = [G\max(G*y,0);zeros(bas_est);zeros(c1_est)];
        myfun = @(Ald) lagrangian_temporal_gradient(Ald,thr^2,y(~mis_data)-b-c1*gd_vec(~mis_data),bas_est,c1_est, c);
        options_dual = optimset('GradObj','On','Display','Off','Algorithm','interior-point','TolX',1e-8);
        ld_in = 10;
        [ld,~,flag] = fmincon(myfun,ld_in,[],[],[],[],0,[],[],options_dual);         %#ok<ASGLU>
        if (flag == -2) || (flag == -3)
            warning('Problem seems unbounded or infeasible. Try a different method.');
        end
        if bas_est; b = c(T+bas_est); end
        if c1_est; c1 = c(end); end
        c = c(1:T);
        sp = G*c;
    case 'cvx'
        onPath = ~isempty(which('cvx_begin'));
        if onPath
            c = zeros(T,1+options.resparse);
            sp = zeros(T,1+options.resparse);
            bas = zeros(1+options.resparse,1);
            cin = zeros(1+options.resparse,1);
            w_ = ones(T,1);
            for rep = 1:options.resparse+1
                [c(:,rep),bas(rep),cin(rep)] = cvx_foopsi(y,b,c1,sn,b_lb,g,w_,~mis_data);
                sp(:,rep) = G*c(:,rep);                
                w_ = 1./(max(sp(:,rep),0) + 1e-8);
            end
            sp(sp<1e-6) = 0;
            c = G\sp;
            b = bas;
            c1 = cin;
        else
            error('CVX does not appear to be on the MATLAB path. It can be downloaded from cvxr.com \n');
        end
    case 'lars'
         Ginv = E*[full(G\speye(T)),ones(T,bas_est),gd_vec*ones(1,c1_est)];
         if bas_est; b = 0; end
         if c1_est; c1 = 0; end    
         [~, ~, spikes, ~, ~] = lars_regression_noise(y(~mis_data)-b_lb*bas_est - b - c1*gd_vec(~mis_data), Ginv, 1, sn^2*(T-sum(mis_data)));
         sp = spikes(1:T);
         b = (spikes(T+bas_est)+b_lb)*bas_est + b*(1-bas_est);
         c1 = spikes(end)*c1_est + c1*(1-c1_est);
         c = G\sp;
    case 'spgl1'
        onPath = ~isempty(which('spgl1'));
        if onPath
            Gx = @(x,mode) G_inv_mat(x,mode,T,g,gd_vec,bas_est,c1_est,E);
            c = zeros(T,1+options.resparse);
            sp = zeros(T,1+options.resparse);
            bas = zeros(1+options.resparse,1);
            cin = zeros(1+options.resparse,1);
            w_ = ones(T,1);
            for rep = 1:options.resparse+1
                if bas_est; b = 0; w_ = [w_;1e-10]; end %#ok<AGROW>
                if c1_est; c1 = 0; w_ = [w_;1e-10]; end %#ok<AGROW>
                options_spgl = spgSetParms('project',@NormL1NN_project ,'primal_norm', @NormL1NN_primal,'dual_norm',@NormL1NN_dual,'verbosity',0,'weights',w_);
                [spikes,r,~,~] = spg_bpdn( Gx, y(~mis_data)-b_lb*bas_est - (1-bas_est)*b-(1-c1_est)*c1*gd_vec(~mis_data), sn*sqrt(T-sum(mis_data)),options_spgl); %#ok<ASGLU>
                c(:,rep) = G\spikes(1:T); %Gx([spikes(1:T);0],1);                                  %% calcium signal
                bas(rep) = b*(1-bas_est) + bas_est*spikes(T+bas_est)+b_lb*bas_est;       %% baseline
                cin(rep) = c1*(1-c1_est) + c1_est*spikes(end);
                sp(:,rep) = spikes(1:T);                                           %% spiking signal
                w_ = 1./(spikes(1:T)+1e-8);
            end
            b = bas;
            c1 = cin;
            %sn = norm(r)/sqrt(T);
        else
            error('SPGL1 does not appear to be on the MATLAB path. It can be downloaded from math.ucdavis.edu/~mpf/spgl1 \n');
        end
end

% constrained.c   = c;
% constrained.b   = b;
% constrained.c1  = c1;
% constrained.g   = g;
% constrained.sn  = sn;
% constrained.sp  = sp;
% constrained.F_est = c;
% constrained.spk = sp;

end

function sn = GetSn(Y,range_ff,method)
    % estimate noise level with a power spectral density method
    L=length(Y);
%         if ~isempty(which('pmtm'))
%             [psd_Y,ff] = pmtm(Y,5/2,1000,1);
%         end
    if ~isempty(which('pwelch'));
        [psd_Y,ff]=pwelch(Y,round(L/8),[],1000,1);
    else
        xdft = fft(Y);
        xdft = xdft(:,1:round(L/2)+1);
        psd_Y = (1/L) * abs(xdft).^2;
        ff = 0:1/L:1/2;
        psd_Y(2:end-1) = 2*psd_Y(2:end-1);
    end
    ind=ff>range_ff(1);
    ind(ff>range_ff(2))=0;
    switch method
        case 'mean'
            sn=sqrt(mean(psd_Y(ind)/2));
        case 'median'
            sn=sqrt(median(psd_Y(ind)/2));
        case 'logmexp'
            sn = sqrt(exp(mean(log(psd_Y(ind)/2))));
    end
end

function g = estimate_time_constants(y,p,sn,lags)
    % estimate time constants from the sample autocovariance function
    lags = lags + p;
    if ~isempty(which('xcov')) %signal processing toolbox
        xc = xcov(y,lags,'biased');
    else
        ynormed = (y - mean(y));
        xc = nan(lags + 1, 1);
        for k = 0:lags
            xc(k + 1) = ynormed(1 + k:end)' * ynormed(1:end - k);
        end
        xc = [flipud(xc(2:end)); xc] / numel(y);
    end
    xc = xc(:);
    A = toeplitz(xc(lags+(1:lags)),xc(lags+(1:p))) - sn^2*eye(lags,p);
    g = pinv(A)*xc(lags+2:end);            
end

function [f,grad] = lagrangian_temporal_gradient(Al,thr,y_raw,bas_flag,c1_flag, c)
    options_qp = optimset('Display','Off','Algorithm','interior-point-convex');
    H = [speye(T),ones(T,bas_flag),gd_vec*ones(1,c1_flag);...
        ones(bas_flag,T),T*ones(bas_flag),(1-gd^T)/(1-gd)*ones(c1_flag,bas_flag);...
        (gd_vec*ones(1,c1_flag))',(1-gd^T)/(1-gd)*ones(c1_flag,bas_flag),(1-gd^(2*T))/(1-gd^2)*ones(c1_flag,c1_flag)];
    Ay = [y_raw;sum(y_raw)*ones(bas_flag);gd_vec'*y_raw*ones(c1_flag)];
    c = quadprog(2*Al(1)*H, [v;zeros(bas_flag+c1_flag,1)]-2*Al(1)*Ay, [-G,sparse(T,bas_flag+c1_flag);sparse(bas_flag+c1_flag,T),-speye(bas_flag+c1_flag)]...
        , [sparse(T,1);-b_lb*ones(bas_flag);zeros(c1_flag)], [], [], [], [], c, options_qp);
    f = v'*c(1:T);    
    grad = [sum((c(1:T)-y_raw + c(T+bas_flag)*bas_flag + c(end)*gd_vec*c1_flag).^2)-thr]; %#ok<NBRAK>
    f = f + Al(:)'*grad;
end

function b = G_inv_mat(x,mode,NT,gs,gd_vec,bas_flag,c1_flag,Emat)
    if mode == 1
        b = filter(1,[1;-gs(:)],x(1:NT)) + bas_flag*x(NT+bas_flag) + c1_flag*gd_vec*x(end);
        b = Emat*b;
       %b = G\x(1:NT) + x(NT+bas_flag)*bas_flag + x(end)*c1_flag;
    elseif mode == 2
        x = Emat'*x;
        b = [flipud(filter(1,[1;-gs(:)],flipud(x)));ones(bas_flag,1)*sum(x);ones(c1_flag,1)*(gd_vec'*x)];
       %b = [G'\x;ones(bas_flag,1)*sum(x);ones(c1_flag,1)*(gd_vec'*x)] ;
    end
end

function [c,b,c1] = cvx_foopsi(y,b,c1,sn,b_lb,g,w,keep)
% implementation of constrained foopsi in CVX
% Written by Eftychios Pnevmatikakis
    if isempty(b)
        bas_est = 1;
    else
        bas_est = 0;
    end
    if isempty(c1)
        c1_est = 1;
    else
        c1_est = 0;
    end
    gd = max(roots([1,-g(:)']));
    T = length(y);
    G = spdiags(ones(T,1)*[-g(end:-1:1)',1],-length(g):0,T,T);
    gd_vec = gd.^((0:T-1)');
    cvx_begin quiet
        variable c2(T)
        if bas_est; variable b; end
        if c1_est; variable c1; end
        minimize(w'*(G*c2))
        subject to
            G*c2>=0;
            norm(y(keep)-c2(keep)-b-c1*gd_vec(keep))<=sqrt(sum(keep))*sn;
            if bas_est; b>=b_lb; end
            if c1_est; c1>=0; end
    cvx_end
    if strcmpi(cvx_status,'Infeasible');
        %disp('Problem is infeasible, adjusting noise value.');
        cvx_begin quiet
            variable c2(T)
            if bas_est; variable b; end
            if c1_est; variable c1; end
            minimize(norm(y(keep)-c2(keep)-b-c1*gd_vec(keep)))
            subject to
                G*c2>=0;
                if bas_est; b>=b_lb; end
                if c1_est; c1>=0; end
        cvx_end
        sn = cvx_optval/sqrt(sum(keep));
    end
    c = c2;
end

function [Ws, lambdas, W_lam, lam, flag] = lars_regression_noise(Y, X, positive, noise)
% run LARS for regression problems with LASSO penalty, with optional positivity constraints
% Author: Eftychios Pnevmatikakis. Adapted code from Ari Pakman
% Input Parameters:
%   Y:           Y(:,t) is the observed data at time t
%   X:           the regresion problem is Y=X*W + noise
%   maxcomps:    maximum number of active components to allow
%   positive:    a flag to enforce positivity
%   noise:       the noise of the observation equation. if it is not
%                provided as an argument, the noise is computed from the
%                variance at the end point of the algorithm. The noise is
%                used in the computation of the Cp criterion.


% Output Parameters:
%   Ws: weights from each iteration
%   lambdas: lambda values at each iteration
%   Cps: C_p estimates
%   last_break:     last_break(m) == n means that the last break with m non-zero weights is at Ws(:,:,n)
verbose=false;
%verbose=true;
k=1;
T = size(Y,2); % # of time steps
N = size(X,2); % # of compartments
maxcomps = N;
W = zeros(N,k);
active_set = zeros(N,k);
visited_set = zeros(N,k);
lambdas = [];
Ws=zeros(size(W,1),size(W,2),maxcomps);  % Just preallocation. Ws may end with more or less than maxcomp columns
r = X'*Y(:);         % N-dim vector
M = -X'*X;            % N x N matrix 
% begin main loop
%fprintf('\n i = ');
i = 1;
flag = 0;
while 1
    if flag == 1;
        W_lam = 0;
        break;
    end
  %  fprintf('%d,',i);
    % calculate new gradient component if necessary
    if i>1 && new && visited_set(new) ==0         
        visited_set(new) =1;    % remember this direction was computed
    end
    % Compute full gradient of Q 
    dQ = r + M*W;
    % Compute new W
    if i == 1
        if positive
            dQa = dQ;
        else
            dQa = abs(dQ);
        end
        [lambda, new] = max(dQa(:));    
        if lambda < 0
            disp('All negative directions!')
            break
        end
    else
        % calculate vector to travel along 
%        disp('calc velocity')
        [avec, gamma_plus, gamma_minus] = calcAvec(new, dQ, W, lambda, active_set, M, positive);        
        % calculate time of travel and next new direction 
        if new==0                               % if we just dropped a direction we don't allow it to emerge 
            if dropped_sign == 1                % with the same sign
                gamma_plus(dropped) = inf;
            else
                gamma_minus(dropped) = inf;
            end 
        end
        gamma_plus(active_set == 1) = inf;       % don't consider active components 
        gamma_plus(gamma_plus <= 0) = inf;       % or components outside the range [0, lambda]
        gamma_plus(gamma_plus> lambda) =inf;
        [gp_min, gp_min_ind] = min(gamma_plus(:));
        if positive
            gm_min = inf;                         % don't consider new directions that would grow negative
        else
            gamma_minus(active_set == 1) = inf;            
            gamma_minus(gamma_minus> lambda) =inf;
            gamma_minus(gamma_minus <= 0) = inf;                
            [gm_min, gm_min_ind] = min(gamma_minus(:));

        end
        [g_min, which] = min([gp_min, gm_min]);
        if g_min == inf;               % if there are no possible new components, try move to the end
            g_min = lambda;            % This happens when all the components are already active or, if positive==1, when there are no new positive directions 
        end
        % LARS check  (is g_min*avec too large?)
        gamma_zero = -W(active_set == 1)  ./ avec;
        gamma_zero_full = zeros(N,k);
        gamma_zero_full(active_set == 1) = gamma_zero;
        gamma_zero_full(gamma_zero_full <= 0) = inf;
        [gz_min, gz_min_ind] = min(gamma_zero_full(:));
        
        if gz_min < g_min            
            if verbose
                fprintf('DROPPING active weight: %d.\n', gz_min_ind)
            end
            active_set(gz_min_ind) = 0;
            dropped = gz_min_ind;
            dropped_sign = sign(W(dropped));
            W(gz_min_ind) = 0;
            avec = avec(gamma_zero ~= gz_min);
            g_min = gz_min;
            new = 0;
            
        elseif g_min < lambda
            if  which == 1
                new = gp_min_ind;
                if verbose
                    fprintf('new positive component: %d.\n', new)
                end
            else
                new = gm_min_ind;
                fprintf('new negative component: %d.\n', new)
            end
        end
                               
        W(active_set == 1) = W(active_set == 1) + g_min * avec;
        
        if positive
            if any (W<0)
                min(W);
                flag = 1;
                %error('negative W component');
            end
        end
        
        lambda = lambda - g_min;
    end
%  Update weights and lambdas 
    lambdas(i) = lambda;
    Ws(:,:,i)=W;
    res = norm(Y-X*W,'fro')^2;
    if lambda ==0 || (new && sum(active_set(:)) == maxcomps) || (res < noise)       
        if verbose
            fprintf('end. \n');
        end
        
        break
    end 
    if new        
        active_set(new) = 1;
    end

    
    i = i + 1;
end% end main loop 
% final calculation of mus
if flag == 0
    if i > 1
        Ws= squeeze(Ws(:,:,1:length(lambdas)));
        w_dir = -(Ws(:,i) - Ws(:,i-1))/(lambdas(i)-lambdas(i-1));
        Aw = X*w_dir;
        y_res = Y - X*(Ws(:,i-1) + w_dir*lambdas(i-1));
        ld = roots([norm(Aw)^2,-2*(Aw'*y_res),y_res'*y_res-noise]);
        lam = ld(intersect(find(ld>lambdas(i)),find(ld<lambdas(i-1))));
        if numel(lam) == 0  || any(lam)<0 || any(~isreal(lam));
            lam = lambdas(i);
        end
        W_lam = Ws(:,i-1) + w_dir*(lambdas(i-1)-lam(1));
    else
        cvx_begin quiet
            variable W_lam(size(X,2));
            minimize(sum(W_lam));
            subject to
                W_lam >= 0;
                norm(Y-X*W_lam)<= sqrt(noise);
        cvx_end
        lam = 10;
    end
else
    W_lam = 0;
    Ws = 0;
    lambdas = 0; 
    lam = 0;
end
end

function [avec, gamma_plus, gamma_minus] = calcAvec(new, dQ, W, lambda, active_set, M, positive)
[r,c] = find(active_set); %#ok<ASGLU>
Mm = -M(r,r);
Mm=(Mm + Mm')/2;
% verify that there is no numerical instability 
eigMm = eig(Mm);
if any(eigMm < 0)
    min(eigMm)
    %error('The matrix Mm has negative eigenvalues')  
    flag = 1;
end
b = sign(W);
if new
    b(new) = sign(dQ(new));
end
b = b(active_set == 1);
avec = Mm\b;
if positive 
    if new 
        in = sum(active_set(1:new));
        if avec(in) <0
            new;
            %error('new component of a is negative')
            flag = 1;
        end
    end
end
one_vec = ones(size(W));
dQa = zeros(size(W));
for j=1:length(r)
    dQa = dQa + avec(j)*M(:, r(j));
end
gamma_plus = (lambda - dQ)./(one_vec + dQa);
gamma_minus = (lambda + dQ)./(one_vec - dQa);
end
