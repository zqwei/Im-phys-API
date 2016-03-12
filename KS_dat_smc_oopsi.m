
function KS_dat_smc_oopsi(nCell)

load('DataListCells.mat')
warning('off', 'all');

dat     = totCell(nCell);
dff     = double(dat.dff);
spk     = double(dat.spk);
caTime  = double(dat.CaTime);
dt      = caTime[1] - caTime[0];
fr      = 1/dt;
T       = length(caTime);
% smc-oopsi
V.fast_iter_max    = 100;
V.smc_iter_max     = 100;
V.dt               = 1/fr;
V.T                = length(dff);
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
P.A     = 50;
P.n     = 1;  
P.k_d   = 200; 

[smc, fast, F, estF] = smc_oopsi(dff, V, P);
smc.F                = F;
smc.F_est            = estF;

[ca_p,peel_p, data]  = peel_oopsi(dff, fr);
peel                 = data;
peel.ca_params       = ca_p;
peel.peel_params     = peel_p;

save(['Smc_oopsi_fit_Cell_' num2str(nCell)], 'smc', 'peel', 'fast')

end
