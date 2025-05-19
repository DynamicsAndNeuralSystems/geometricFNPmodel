%% Set topology input parameters

topology = struct;

topology.T = 0.07;
topology.L = 0.4;

%% Set homogeneous input parameters

homparam = struct;

homparam.gamma = 116;
homparam.r = 0.086;
homparam.nu0 = 0.756;

%% Set topology grid size and timestep in accordance to homogeneous parameters

topology.Nx = 200;
dx = topology.L / topology.Nx;
dt = dx / (homparam.r * homparam.gamma * sqrt(2) * 2);
topology.Nt = ceil(topology.T / dt);
clear dx dt;

%% Set homogeneous and heterogeneous input parameters

hetparam_hom = struct;
hetparam_hom.m = 0;

hetparam_het = struct;
hetparam_het.m = 1;
hetparam_het.c = (homparam.r)^2;
hetparam_het.tau = 0;
hetparam_het.a = [0.15; 0.15];
hetparam_het.b = [0.25; 0.25];
hetparam_het.sigmaeps = 0.002;

%% Set stimulation parameters

stim = struct;

stim.stimnum = 1;
stim.sigma = [0.004, 0.0006];

stimt = zeros(1,1);
stimt(1) = 0.02;
stimI = zeros(1,1);
stimI(1) = 1;
stimR = zeros(2, 1);
stimR(:, 1) = [0.2, 0.2];

stim.stimR = stimR;
stim.stimt = stimt;
stim.stimI = stimI;


clear stimt stimI stimR;

