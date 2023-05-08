% This script generates the input files calculating linear abs for NISE
% It will also run the NISE calculations provided that the correct path to
% NISE program is ndicated.

% Input:
% - Site energies (text)
% - Protein structure
% Output:
% - energy trajectories
% - dipole trajectories
% - NISE input files (translation, Absorption, 2DES, ...)

%% Files
f_pdb = {fullfile('pdb','Chlab.pdb'),'Y'}; % Structure pdb
f_site = fullfile('Energy','Chlab.txt'); % Site energies

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

%% Parameters
sigma = [60 140]; % Disorders (static & dynamic) [cm-1]
tauc = [150 1000000]; % Correlation time [fs]
dt = 1; %Time step for trajectories [fs]
Nstep = 20000; % Number of time steps
taudeph = 150; % Dephasing time [fs]
Tw = 0; % Waiting time for 2DES [fs]
T = 300; % Temperature (K)

E0 = load(f_site); % Site energies [cm-1]
N = length(E0); % Number of chromophores
maxfreq = max(E0) + 1000;
minfreq = min(E0) - 1000;

%% NISE input parameters
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
nise1D.RunTimes = [1.5*round(taudeph/dt) 0 1.5*round(taudeph/dt)];
nise1D.Singles = N;

% NISE input for 2D
nise2D = nise1D;
nise2D.Technique = '2DUVvis';
nise2D.RunTimes = [1.5*round(taudeph/dt) Tw 1.5*round(taudeph/dt)];

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

%% Generate energy trajectories
dE = odam_trajectory(E0,0:dt:Nstep*dt,sigma,1./tauc); % Energy fluctuation [cm-1]
E = E0 + dE;

%% Calculate coupling
atom = import_pdb(f_pdb{1},f_pdb{2});
C = calc_coupling(atom);
mu = C.Dvec.*sqrt(C.D);
box = [max([atom.x])-min([atom.x]), ...
    max([atom.y])-min([atom.y]), ...
    max([atom.z])-min([atom.z])];
boxmax = ceil(max(box)); % Cubic box size containing all atoms

%% Write to files
% Generate energy and dipole
fid_ham = fopen([f_ham '.bin'],'w');
fid_dp = fopen([f_dp '.bin'],'w');
fid_pos = fopen([f_pos '.bin'],'w');
for nt = 1:Nstep
    H = tril(diag(E(:,nt))+C.V);
    H = H(:);
    H(H==0) = []; % Convert to one line
    %     fprintf(fid_ham,'%d %f',nt,H); % write to Hamiltonian file
    %     fprintf(fid_ham,'\n');
    %     fprintf(fid_dp,'%d %f',nt, mu); % write to dipole file
    %     fprintf(fid_dp,'\n');
    %     fprintf(fid_pos,'%d %f',boxmax, C.Center); % write to position file
    %     fprintf(fid_pos,'\n');
    fwrite(fid_ham, [nt; H],'float32');
    fwrite(fid_dp, [nt; mu(:)],'float32');
    fwrite(fid_pos,[boxmax; C.Center(:)],'float32');
end
fclose(fid_ham);
fprintf('Hamiltonian file generated\n');
fclose(fid_dp);
fprintf('Dipole file generated\n');
fclose(fid_pos);
fprintf('Position file generated\n');

% Generate NISE input files
% TXT to BIN Translation
writeinput(niseTra,f_iTra);

% Absorption
writeinput(nise1D,f_i1D);

% 2D
writeinput(nise2D,f_i2D);

% Dif
writeinput(niseDif,f_iDif);

% Pop
writeinput(nisePop,f_iPop);

% Analyse
writeinput(niseAnalyse,f_iAnalyse);

% LD
writeinput(niseLD,f_iLD);

% CD
writeinput(niseCD,f_iCD);

% Luminesence
writeinput(niseLum,f_iLum);

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
