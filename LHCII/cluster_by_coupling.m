function [G,ig] = cluster_by_coupling(V,Vc,E,Ec)
% CLUSTER_BY_COUPLING - Partition the Hamiltonian into hard clusters
% based on the coupling between them

% Updated 20230414: Now loop until Vs reach cutoff, despite all sites are
%                   assigned

% Syntax
%   G = cluster_by_coupling(V,Vc)
%   G = cluster_by_coupling(V,Vc,E,Ec)
%
% Input
%    V - couplings between sites (N-by-N matrix)
%    Vc - cutoff value for the coupling
%
%    (Optional)
%    E - site energies (N-by-1)
%    Ec - maximal site energy difference to include in the same domain
%
% Output
%    G - cluster index for each site (N-by-1)

N = length(V);

% Get the lower triangle of V (V should be symmetric)
Vt = abs(tril(V));

% Remove couplings between sites with large energy difference (dE > Ec)
if exist('Ec','var')
    Em = repmat(E,1,N);
    En = repmat(E',N,1);
    dE = Em-En;
    Vt(abs(dE) > Ec) = 0;
end

% Sort the couplings from largest to smallest
% Vs is an N*N-element vector of sorted couplings Vmn
% ix is an 2-column matrix showing the coupled pairs
[Vs, ix] = sort(Vt(:),'desc');
[ia, ib] = ind2sub([N,N],ix);
ix = [ia, ib];

% G is the cluster index map
G = zeros(N,1);

% Assign the same cluster to each coupled pair starting with the strongest
% Loop until cutoff Vc is reached
k = 1; % loop index
g = 0; % cluster counter
while Vs(k) > Vc
    ik = ix(k,:); % get the k-th pair
    ig = G(ik)>0; % read the assigned clusters for pair k
    if all(ig)
        % clusters are assigned to each site - combine clusters
        g1 = min(G(ik)); g2 = max(G(ik));
        G(G==g2) = g1;
    elseif any(ig)
        % cluster is assigned to one site - assign the same to the other        
        G(ik(~ig)) = G(ik(ig));
    else
        % no cluster is assigned yet - make a new cluster
        g = g + 1;
        G(ik) = g;
        % expand cluster
        ij = Vt(:,ik(1)) > Vc & (abs(E-E(ik(1))) < Ec);
        G(ij) = g;
        ij = Vt(:,ik(2)) > Vc & (abs(E-E(ik(2))) < Ec);
        G(ij) = g;
    end
    k = k + 1;
end

% After exiting the loop, if there are remaining sites, assign
% each of them to its own single-site cluster
orphans = G==0;
G(orphans) = g+1:g+numel(find(orphans));

% Replace cluster IDs to sequential numbers
[GID,~,ic] = unique(G);
SEQ = 1:numel(GID);
G = SEQ(ic)';

% Rearrange Hamiltonian
Clusters = unique(G)';
ig = false(N,numel(Clusters));
for k = Clusters
    ki = G==k;     % sites in cluster k
    ig(:,k) = ki;
end

% Make G starts at zero
seg_name = unique(G,'stable');
seg_remap = 0:(length(seg_name)-1);
G2 = G;
for i = 1:length(G)
    G2(i) = seg_remap(seg_name==G(i));
end
G = G2;
end
