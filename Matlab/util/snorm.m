function [snorm, matrixnorm, C] = snorm(Lambda, s)
N = sqrt(size(Lambda,1)/2)-1;

syms x
leg_zeros = double(vpasolve(legendreP(N+1,x) ==0));
thph = quad_points_unit_sphere(leg_zeros,'spherical');
th = thph(:,1);
phi = thph(:,2);

%% gauss-legendre weights
[~,alpha]=gauss_legendre(N+1);
alpha = pi/(N+1)*alpha';

%% spherical harmonics
% initialize matrix of spherical harmonics
Y = zeros(2*(N+1)^2,(N+1)^2); 

% Compute Y for each degree and order evaluated in quadrature points
for n = 0:N
    for m = -n:n
        index = n^2+m+(n+1);
        Y(:, index) = harmonicY(n,m,th,phi,'norm',true,'phase',false);
    end
end

%% Compute BY_n^m in each quadrature point
BY = zeros(2*(N+1)^2,(N+1)^2); 
for i = 1:(N+1)^2
    BY(:, i) = Lambda*Y(:,i);
end

%% s-norm weights
w2 = @(n,s) (1+n^2)^s;

%% Compute sum
C = zeros((N+1)^2,(N+1)^2);
alphakl = repmat(alpha,[2*(N+1),1]); 
for n =0:N % row
    for m =-n:n % row
        for nn = 0:N % column
            for mm = -nn:nn  % column    
                index = n^2+m+n+1; % rows in matrix
                iindex = nn^2+mm+nn+1; % columns in matrix
                mneg_index = n^2-m+(n+1);
                C(index,iindex) = w2(nn,-s)*sum(alphakl.*BY(:,iindex).*Y(:,mneg_index));
            end
        end
    end
    n
end

%% 
snorm = norm(C);
matrixnorm = norm(Lambda);
