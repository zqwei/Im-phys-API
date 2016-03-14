warning('off', 'all');
% set simulation metadata
T       = 1000; % # of time steps
V.dt    = 1/8;  % time step size
% initialize params
P.a     = 1;    % observation scale
P.b     = 0;    % observation bias
tau     = 1.5;    % decay time constant
P.gam   = 1-V.dt/tau; % C(t) = gam*C(t-1)
P.lam   = 0.1;  % firing rate = lam*dt
P.sig   = 0.1;  % standard deviation of observation noise 
% simulate data
N = poissrnd(P.lam*V.dt*ones(T,1)); % simulate spike train
C = filter(1,[1 -P.gam],N);         % calcium concentration
F = P.a*C+P.b + P.sig*randn(T,1);   % observations
% smc-oopsi
fr                 = 8; % Hz
V.fast_iter_max    = 1;
V.smc_iter_max     = 100;
V.dt               = 1/fr;
V.T                = length(F);
V.preprocess       = 1;
V.n_max            = 2;
V.Ncells           = 1;
V.Npixels          = 1;
V.T_o              = V.T;
V.x                = ones(1,V.T);
V.est_c            = 1;
V.est_t            = 1;
V.est_n            = 1;
V.est_h            = 0;
V.est_F            = 1;
V.scan             = 0;
V.Nparticles       = 99;
V.Nspikehist       = 0;
V.condsamp         = 1;
V.ignorelik        = 1;
V.use_true_n       = 0;
V.freq             = 1;
Pest.A     = 50;
Pest.n     = 1;  
Pest.k_d   = 200; 

tvec=0:V.dt:(T-1)*V.dt;

% figure;
% [~, ~, data] = peel_oopsi(F', fr);
% subplot(211)
% hold on
% plot(tvec,F)
% plot(data.tim,data.model)
% title('peel')
% subplot(212)
% hold on
% stem(tvec,N); 
% stem(data.tim,data.spiketrain)
% title('peel')

figure;
[~, ~, data] = peel_nl_oopsi(F', fr);
subplot(211)
hold on
plot(tvec,F)
plot(data.tim,data.model)
title('peel')
subplot(212)
hold on
stem(tvec,N); 
stem(data.tim,data.spiketrain)
title('peel')

