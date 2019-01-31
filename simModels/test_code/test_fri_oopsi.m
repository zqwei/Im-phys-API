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

tvec=0:V.dt:(T-1)*V.dt;
fr   = 1/V.dt;
figure;
fri = fri_oopsi(F,tau, fr);
subplot(211)
hold on
plot(tvec,F, 'linewid', 1)
plot(tvec,fri.F_est, 'linewid', 1)
% title('fri')
xlim([0 100])
ylabel([0 1])
xlabel('time (sec)')
ylabel('normalized F')
box off
legend('raw', 'fri')


subplot(212)
hold on
stem(tvec,N); 
stem(fri.spk, ones(length(fri.spk),1))
% title('fri')
xlim([0 100])
ylabel([0 1])
xlabel('time (sec)')
ylabel('normalized # spk')
box off
setPrint(8, 6*2, 'fri_oopsi', 'png')
setPrint(8, 6*2, 'fri_oopsi')



