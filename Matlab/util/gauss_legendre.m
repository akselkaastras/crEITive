function [x w]=gauss_legendre(n)
% [x w] = GAUSS_LEGENDRE(n) computes the nodes and weights for the
% Gauss-Legendre quadrature rule of order n on [-1,1] with the Golub-Welsch
% algorithm. That is, writing the recursion formula on the Legendre monic
% polynomials Q_n:
% For m>=1, Q_{m+1} + (B_m - X)Q_m + A_m Q_{m-1} = 0
%
% The nodes are given by the eigenvalues of the Jacobi tridigonal matrix:
%
% (   B_0     sqrt(A_1)     0          ...             0       )
% ( sqrt(A_1)   B_1      sqrt(A_2)     ...             0       )
% (    0      sqrt(A_2)    B_3         ...             0       )
% (   ...       ...        ...         ...       sqrt(A_{n-1}) )
% (    0         0         ...     sqrt(A_{n-1})    B_{n-1})   )
%
% Considering v_p, 1<=p<=n, normalized eigenvectors of the matrix, the
% weights are given by:
% mu_0 v_p(1)^2
% Where mu_0 is the integral of the weight function omega for the
% orthogonality of the polynomials.
%
% In our particular case of Legendre polynomials:
% (m+1) P_{m+1} - (2m+1) X P_m + m P_{m-1} = 0
% We have :
% A_m = m^2/(4m^2-1)
% B_m = 0
% omega = 1
% mu_0 = 2
%
% x is the row vector of nodes.
% w is the row vector of weights.

% Jacobi matrix
An=(1:n-1).';
An=An./sqrt((4*An.^2-1));
J=spdiags([[An;0],[0;An]],[-1,1],n,n);

% Gauss-Legendre nodes and weights
[w x]=eig(full(J));
clear J;
x=diag(x).';
w=2*w(1,:).^2;

% Symmetrize
x=(x-fliplr(x))/2;
w=(w+fliplr(w))/2;