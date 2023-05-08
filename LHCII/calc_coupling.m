function out = calc_coupling(atom)
% Calculate coupling between Chl molecules
%
% Input parameters:
%  atom - struct containing atom coordinates (MG, NB, ND)
%
% Output parameters (in struct out):
%  N - number of molecules
%  molid - molecule IDs
%  R - distances
%  V - coupling energies
%  k2 - orientation factors
%  D - dipole strength for each Chl
%  Dvec - normalized dipole vector
%  Center - center of the pophyrin ring

out = struct;

% D - dipole strength (Debye^2)
% Muh et al 2010 JPC
DS.CHL = 15; % Chl b
DS.CLA = 21; % Chl a
DS.PHO = 13; % Pheophytin (Raszewski 2008)
DS.BCL = 44; % BChla
%  n - refractive index
n = 1.4;

% Influence of the medium
C = 1/n^2;
% C = (n^2+2)^2/(9*n^2);

% Split MG, NB and ND atom tables
atom_mg = atom(strcmp({atom.type},'MG'));
atom_na = atom(strcmp({atom.type},'NA'));
atom_nb = atom(strcmp({atom.type},'NB'));
atom_nc = atom(strcmp({atom.type},'NC'));
atom_nd = atom(strcmp({atom.type},'ND'));

% Molecule IDs
molid = [atom_na.molid]';
out.molid = molid;
out.chain = [atom_na.chain]';
out.resname = {atom_na.resname}';

% Initialize matrices
N = length(molid); % number of pigments
out.N = N;
out.R = zeros(N,N);  % distances between pigments (nm)
out.V = zeros(N,N);  % coupling energies (cm^-1)
out.k2 = zeros(N,N); % orientation factors
out.Dvec = zeros(N,3); % x y z components of dipole (normalized) [Angstrom]
%% Dipole strength
D = zeros(N,1);
for k = 1:N
    D(k) = DS.(atom_na(k).resname);
%     if strcmp(atom_na(k).chain,'A') || strcmp(atom_na(k).chain,'B') || strcmp(atom_na(k).chain,'C') || strcmp(atom_na(k).chain,'D')
%         D(k) = D(k)*2; %dipole in core is 7/3 times in LHCII
%     end
end
out.D = D;

%% Calculate transition dipole moments
mu = zeros(N,3);

for k = 1:N
    mu(k,:) = [atom_nd(k).x atom_nd(k).y atom_nd(k).z] - ...
        [atom_nb(k).x atom_nb(k).y atom_nb(k).z];
end

% out.mu = mu;

%% Loop over molecules
for a = 1:N
        % Centre
%     C_a1 = [atom_mg(a).x atom_mg(a).y atom_mg(a).z];
    Center(a,:) = mean([atom_na(a).x atom_na(a).y atom_na(a).z;
                atom_nb(a).x atom_nb(a).y atom_nb(a).z;
                atom_nc(a).x atom_nc(a).y atom_nc(a).z;
                atom_nd(a).x atom_nd(a).y atom_nd(a).z]);
end
out.Center = Center;
for a = 1:N
    % Centre
    C_a = Center(a,:);
    
    % Transition dipole moment
    mu_a = mu(a,:);
    
    % normalize (unit vectors)
    norm_a = norm(mu_a);
    mn_a = mu_a ./ norm_a;
    out.Dvec(a,:) = mn_a;
    
    for b = 1:N
        if b~=a
            % Centre
            C_b = Center(b,:);
            
            % Transition dipole moment
            mu_b = mu(b,:);
            
            % normalize (unit vectors)
            norm_b = norm(mu_b);
            mn_b = mu_b ./ norm_b;
            
            % Dot product
            ab = (dot(mn_a,mn_b));
            
            % Center-to-center vector r
            r = (C_b - C_a)./10;
            R = norm(r); % Distance in nm
            rn = r ./ R; %            
            out.R(b,a) = R;
            
            % Dot products of r with the transition dipole moments
            ra = (dot(rn,mn_a));
            rb = (dot(rn,mn_b));
            
            % Orientation factor
            kappa = ab - 3*ra*rb;
            out.k2(b,a) = kappa^2;
            
            % Coupling strength 
%             Dab = sqrt(D(a)*D(b)); 
            % (point-dipole approximation)
%             out.V(b,a) = 5.04*C*Dab*kappa/R^3; % cm^-1

            % Dipole approximation (Shibata et al 2010 JPC)
            mu_a = mn_a*sqrt(D(a));
            mu_b = mn_b*sqrt(D(b));
            out.V(b,a) = 5.04*C*(dot(mu_a,mu_b)/R^3 ...
                         - 3*(dot(mu_a,r)*dot(mu_b,r))/R^5);           
        end
    end
end

end
