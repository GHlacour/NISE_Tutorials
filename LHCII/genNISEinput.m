% This script generates the input files calculating linear abs for NISE
% It will also run the NISE calculations provided that the correct path to
% NISE program is indicated.

% Input:
% - Site energies (text)
% - Protein structure
% - Several parameters
% Output:
% - energy trajectories
% - dipole trajectories
% - NISE input files (translation, Absorption, 2DES, ...)

%% Files and parameters
f_pdb = {fullfile('pdb','5xnm.pdb'),'Y'}; % Structure pdb and chain name
f_site = fullfile('Energy','LHCIImon.txt'); % Site energies

f_ham = 'Energy'; % Hamiltonian input for NISE
f_dp = 'Dipole'; % Dipole input for NISE
f_pos = 'Position'; % Position input for NISE
f_i1D = 'input1D'; % NISE absorption calculation parameters
f_iTra = 'inpTra'; % NISE translate parameters
f_i2D = 'input2D'; % NISE 2D calculation parameters
f_iCG2D = 'inputCG2D'; % NISE coarse grained 2D calculation
f_iDif = 'inputDif'; % NISE Diffusion
f_iPop = 'inputPop'; % NISE population transfer
f_iAnalyse = 'inputAnalyse'; % NISE Analyze
f_iLum = 'inputLum'; % NISE Luminescence
f_iCD = 'inputCD'; % NISE circular dichroism
f_iLD = 'inputLD'; % NISE linear dichroism
f_iDOS = 'inputDOS'; % NISE density of states
f_MCFRET = 'inputMCFRET'; % NISE MCFRET rates

% Parameters
sigma = [180 60] % Disorders (dynamic & static) [cm-1]
tauc = [150 10000] % Correlation time [fs]
dt = 1 %Time step for trajectories [fs]
Nstep = 300000 % Number of time steps
taudeph = 75 % Pure Dephasing time [fs]
Tw = 0 % Waiting time for 2DES [fs]
T = 300 % Temperature (K)
tmax = 3*min([tauc taudeph]) % t1 and t3 max [fs]
freq_margin = 1000 % Frequency axis = [minfreq-margin maxfreq+margin] [cm-1]
cluster_coup_cut = 20 % Coupling cutoff for clustering [cm-1]
cluster_freq_cut = 300 % Frequency cutoff for clustering [cm-1]

%% Initializing variables
E0 = load(f_site); % Site energies [cm-1]
fprintf('Energy loaded from %s\n',f_site);
N = length(E0); % Number of chromophores
maxfreq = max(E0) + freq_margin;
minfreq = min(E0) - freq_margin;
t1 = 0:dt:(Nstep-1)*dt; % Time axis for trajectories
sigma = repmat(sigma,N,1);

%% Generate energy trajectories
dE = odam_trajectory(E0,t1,sigma,1./tauc); % Energy fluctuation [cm-1]
E = E0 + dE;
fprintf('Trajectory calculated\n');

%% Calculate coupling
atom = import_pdb(f_pdb{1},f_pdb{2});
C = calc_coupling(atom);
mu = C.Dvec.*sqrt(C.D);
box = [max([atom.x])-min([atom.x]), ...
    max([atom.y])-min([atom.y]), ...
    max([atom.z])-min([atom.z])];
boxmax = ceil(max(box)); % Cubic box size containing all atoms
fprintf('Coupling calculated\n');
cluster_index = cluster_by_coupling(C.V,cluster_coup_cut,E0,cluster_freq_cut)';

%% Write to files

fprintf('Write to files...\n');
writematrix(cluster_index,'cluster.txt');
% Generate energy and dipole
fid_ham = fopen([f_ham '.bin'],'w');
fid_dp = fopen([f_dp '.bin'],'w');
fid_pos = fopen([f_pos '.bin'],'w');

Hmask = true(N);
Hmask = tril(Hmask); % Create a filter to get lower triagular elements
for nt = 1:Nstep
    H = diag(E(:,nt))+C.V;
    H = H(Hmask);
    fwrite(fid_ham, [nt; H],'float32');
    fwrite(fid_dp, [nt; mu(:)],'float32');
    fwrite(fid_pos,[boxmax; C.Center(:)],'float32');
end

fclose(fid_ham);
format shortG;
disp('Hamiltonian file generated');
% disp('Last Hamiltonian:');
% disp(diag(E(:,end))+C.V);
% disp('Triangular Hamiltonian:');
% disp(H');
fclose(fid_dp);
disp('Dipole file generated');
fclose(fid_pos);
disp('Position file generated');

%% Generate NISE input files
% For translate between bin and txt
niseTra.InputEnergy = [f_ham '.bin'];
niseTra.InputDipole = [f_dp '.bin'];
niseTra.OutputEnergy = [f_ham '.txt'];
niseTra.OutputDipole = [f_dp '.txt'];
niseTra.Singles = N;
niseTra.Doubles = 0;
niseTra.Skip = 'Doubles';
niseTra.Length = Nstep;
niseTra.InputFormat = 'GROBIN';
niseTra.OutputFormat = 'GROASC';

% NISE input parameters for absorption
nise1D.Propagation = 'Coupling';
nise1D.Couplingcut = 0;
nise1D.Threshold = 0.001;
nise1D.Hamiltonianfile = [f_ham '.bin'];
nise1D.Dipolefile = [f_dp '.bin'];
nise1D.Length = Nstep;
nise1D.Samplerate = 20;
nise1D.Lifetime = taudeph;
nise1D.Timestep = dt;
nise1D.Trotter = 1;
nise1D.Anharmonicity = 100;
nise1D.Format = 'Dislin';
nise1D.MinFrequencies = [minfreq minfreq minfreq];
nise1D.MaxFrequencies = [maxfreq maxfreq maxfreq];
nise1D.Technique = 'Absorption';
nise1D.FFT = 2048;
nise1D.RunTimes = [round(tmax/dt) 0 round(tmax/dt)];
nise1D.Singles = N;

% NISE input for 2D
nise2D = nise1D;
nise2D.Technique = '2DUVvis';
nise2D.RunTimes = [round(tmax/dt) Tw/dt round(tmax/dt)];

% NISE input for CG-2D
niseCG2D = nise2D;
nise.Technique = 'CG_2DES';

% NISE input for diffusion
niseDif = nise1D;
niseDif.Technique = 'Dif';
niseDif.Positionfile = [f_pos '.bin'];

% NISE input for population transfer
nisePop = nise1D;
nisePop.Technique = 'Pop';
nisePop.Propagation = 'Sparse';

% NISE input for analyze
niseAnalyse = nise1D;
niseAnalyse.Technique = 'Analyse';

% NISE input for luminescence
niseLum = nise1D;
niseLum.Technique = 'Luminescence';
niseLum.Temperature = T;

% NISE input for CD
niseCD = nise1D;
niseCD.Technique = 'CD';
niseCD.Positionfile = [f_pos '.bin'];

% NISE input for LD
niseLD = nise1D;
niseLD.Technique = 'LD';

% NISE input for DOS
niseDOS = nise1D;
niseDOS.Technique = 'DOS';

% NSIE input for MCFRET
niseMC = nise1D;
niseMC.Technique = 'MCFRET';
niseMC.Temperature = T;
niseMC.Project = '';
niseMC.Sites = sprintf('%d\n',[N cluster_index]);
% niseMC.Sites = sprintf('%d\n',[N 0 1]);

writeinput(niseTra,f_iTra);
writeinput(nise1D,f_i1D);
writeinput(nise2D,f_i2D);
writeinput(niseDif,f_iDif);
writeinput(nisePop,f_iPop);
writeinput(niseAnalyse,f_iAnalyse);
writeinput(niseLD,f_iLD);
writeinput(niseCD,f_iCD);
writeinput(niseLum,f_iLum);
writeinput(niseMC,f_MCFRET);

%% If NISE is installed, run the calculations
niseDir = '~/NISE_2017/bin/';
if exist(niseDir,'dir') && isunix
    runornot = questdlg('Run NISE calculation? It may take a while.');
    system(sprintf('%s%s %s',niseDir,'translate',f_iTra));
    system(sprintf('%s%s %s',niseDir,'NISE',f_i1D));
    system(sprintf('%s%s %s',niseDir,'NISE',f_iLum));
    system(sprintf('%s%s %s',niseDir,'NISE',f_iCD));
    system(sprintf('%s%s %s',niseDir,'NISE',f_iLD));
    system(sprintf('%s%s %s',niseDir,'NISE',f_i2D));
    system(sprintf('%s%s %s',niseDir,'2DFFT',f_i2D));
end
