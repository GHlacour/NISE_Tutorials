K = [-0.1 0.01;
     0.1  -0.011];
[v,d] = eig(K);
iv = inv(v);

P0 = [1 0;
      0 1];
t = 0:1:99;
P = zeros([length(t),numel(P0)]);
for i = 1:length(t)
    x = v*diag(exp(diag(d)*t(i)))/v*P0;
%     x = expm(K*t(i))*P0;
    P(i,:) = x(:);
end
figure; plot(P)

dat = load('KPop.dat');
figure; plot(dat(:,2:end));