function [DNnoisy,NOISE,noiselevelrel]=add_noise_dnmap(DN,noiselevelabs)
% DNnoisy=ADD_NOISE_DNMAP(DN) add a standard normal gaussian noise to the
% Dirichlet-to-Neumann map matrix with a noise level equal to noiselevel.

DNc=whos('DN');
DNc=DNc.complex;

NOISE=randn(size(DN)); % Gaussian

%%%%%%
% Make sure NOISE is in Y (has the same properties as Dirichlet to Neumann
% maps
%%%%%%
%% NOISE(1) = 0

for i = 1:size(DN,1)
    NOISE(i,:) = NOISE(i,:)-sum(NOISE(i,:))/size(DN,1);
end
for i = 1:size(DN,1)
    NOISE(:,i) = NOISE(:,i)-sum(NOISE(:,i))/size(DN,1);
end

%% int_{\partial /Omega} NOISE f = 0
N = sqrt(size(DN,1)/2)-1;
[~,alpha]=gauss_legendre(N+1);
alpha = pi/(N+1)*alpha;
c = repmat(alpha,[1,2*(N+1)]);
T = sum(c);
C = repmat(alpha,[2*(N+1)^2,2*(N+1)]);

NOISE = (NOISE-1/T*C*NOISE);

%% Add noise 
delta = noiselevelabs/snorm(NOISE,1/2);
DNnoisy = DN+delta*NOISE;
% Compute relative noise level
noiselevelrel = noiselevelabs/(snorm(DN,1/2));
