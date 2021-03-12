function [xix,xiy,xiz,qhat] = read_qhat(filename,opt)
% [xix,xiy,xiz,qhat] = READ_QHAT(filename)
% [xix,xiy,xiz,qhat] = READ_QHAT(filename,'ndgrid')
% [xix,xiy,xiz,qhat] = READ_QHAT(filename,'meshgrid')
%
% Reads qhat computed by the C++ codes eit, fast_eit, eit_cmplx,
% fast_eit_cmplx.
%
% A qhat file is composed of a header given by a positive integer n
% followed on the same line by a real number xim. It means that qhat is
% computed in the cube [-xim,xim]^3 with n regularly spaced discretization
% points of [-xim,xim].
%
% The following of the file is composed of five columns. The first three
% columns give the cartesian coordinates of the points where qhat is
% computed and the last two columns give the values of the real and
% imaginary part of qhat.
%
% The p-th line of the qhat file, where p = n^2*(i-1) + n*(j-1) + k + 1 is
% composed of:
% xix(i) xiy(j) xiz(k) qhat(xix(i),xiy(j),xiz(k)) (real and imaginary part)
%
% If no option is specified or if 'ndgrid' option is given, [xix,xiy,xiz]
% will be given by ndgrid(linspace(-xim,xim,n)) and qhat will be stored in
% a 3 dimensional array giving the values of qhat at the ndgrid points.
%
% If 'meshgrid' option is given, [xix,xiy,xiz] will be given by
% meshgrid(linspace(-xim,xim,n)) and qhat will be stored in a 3 dimensional
% array giving the values of qhat at the meshgrid points.

% Option
if nargin==1
    opt='ndgrid';
elseif nargin==2
    if ~(strcmp(opt,'ndgrid')||strcmp(opt,'meshgrid'))
        error('Bad option, choose nothing, ''ndgrid'' or ''meshgrid''');
    end
else
    error('Bad number of input arguments');
end

% Open file
fid=fopen(filename,'r');

% Check error
if fid==-1
    error(['Error when opening file ',filename]);
end

% Read the first line
n=fscanf(fid,'%u',1);fgetl(fid);

% Read the following of the file
qhat=fscanf(fid,'%e',[5,n^3]).';

% Close file
fclose(fid);

% Post processing
xix=reshape(qhat(:,1),[n,n,n]);
xiy=reshape(qhat(:,2),[n,n,n]);
xiz=reshape(qhat(:,3),[n,n,n]);
qhat=reshape(complex(qhat(:,4),qhat(:,5)),[n,n,n]);

switch opt
    case 'ndgrid'
        xix=permute(xix,[3,2,1]);
        xiy=permute(xiy,[3,2,1]);
        xiz=permute(xiz,[3,2,1]);
        qhat=permute(qhat,[3,2,1]);
    case 'meshgrid'
        xix=permute(xix,[2,3,1]);
        xiy=permute(xiy,[2,3,1]);
        xiz=permute(xiz,[2,3,1]);
        qhat=permute(qhat,[2,3,1]);
    otherwise
        error('Bad option, choose nothing, ''ndgrid'' or ''meshgrid''');
end