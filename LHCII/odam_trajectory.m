function dE = odam_trajectory(E,t,sigma,Lambda)
% This function results in fluctuation of the energies E with time t,
% based on the overdamped brownian oscillator model
% E, sigma and Lambda should be in columns
% In the case one energy fluctuates with multiple modes, the modes are in
% the same row

% Check errors
if ~iscolumn(E)
    error('Energy should be column vector.\n');
end
if size(sigma,1)~=1 && size(sigma,1)~=length(E)
    error('Disorders for one mode should be in a column.\n If there are multiple modes for a site, put them in the same row.\n');
end
if size(Lambda,1)~=1 && size(Lambda,1)~=length(E)
    error('Correlation time for one mode should be in a column.\n If there are multiple modes for a site, put them in the same row.\n');
end
if size(sigma,2)~=size(Lambda,2)
    error('Make sure sigma and lambda have the same dimensions.\n');
end
if ~isrow(t)
    error('Time should be a row vector.\n');
end

% Calculate
dE = zeros(length(E),length(t));
dE(:,1) = sum(sigma.*randn(length(E),size(sigma,2)),2);
% if any(Lambda==0)
%     sigma(Lambda==0) = [];
%     Lambda(Lambda==0) = [];
% end
for nt = 2:length(t)
    dE(:,nt) = dE(:,nt-1).*prod(exp(-Lambda*(t(nt)-t(nt-1))),2) + sum(sigma.*randn(length(E),size(sigma,2)).*sqrt(1-exp(-2*Lambda*(t(nt)-t(nt-1)))),2);
end
end