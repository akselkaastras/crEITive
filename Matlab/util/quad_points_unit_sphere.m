function varargout = quad_points_unit_sphere(t,coord)
% [theta,phi] = QUAD_POINTS_UNIT_SPHERE(t)
%           x = QUAD_POINTS_UNIT_SPHERE(t,coord)
%
% Computes the quadrature points on the sphere, given in spherical
% coordinates by:
% x_{kl}=(theta_k,phi_l)
% That is in cartesian coordinates by:
% x_{kl}=(sin(theta_k)cos(phi_l),sin(theta_k)sin(phi_l),cos(theta_k))
% Where:
% theta_k=arccos(t_k) and phi_l=l*pi/n for k=1,...,n and l=0,...,2n-1
% t_k being the zeros of the Legendre polynomial P_n
%
% t is a vector.
%
% The coord option is either 'cartesian' or 'spherical'. If not specified,
% the the function gives only the angles theta and phi of the quadrature
% points.
%
% If 'cartesian' or 'spherical' is chosen, the quadrature points are
% columnwise organised by increasing first the integer k and then the
% integer l, that is, the column given by a Matlab reshape of the matrix
% X=x_{kl}
%
% x_{1,0}
% x_{2,0}
% ...
% x_{n,0}
% x_{1,1}
% x_{2,1}
% ...
% x_{n,1}
% ...
% ...
% ...
% x_{1,2n-1}
% x_{2,2n-1}
% ...
% x_{n,2n-1}

% n
n=length(t);

% Polar angles
theta=acos(t(:));

% Azimuthal angles
phi=pi/n*(0:2*n-1);

if nargin==1
    varargout{1}=theta;
    varargout{2}=phi.';
elseif nargin==2
    switch coord
        case 'cartesian'
            varargout{1}=[reshape(sin(theta)*cos(phi),2*n^2,1),...
                reshape(sin(theta)*sin(phi),2*n^2,1),...
                repmat(cos(theta),2*n,1)];
        case 'spherical'
            varargout{1}=[repmat(theta,2*n,1),reshape(repmat(phi,n,1),2*n^2,1)];
        otherwise
            error('Bad option, choose nothing, cartesian or spherical');
    end
else
    error('Bad number of input arguments');
end