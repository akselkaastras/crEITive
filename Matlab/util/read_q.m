function varargout = read_q(filename)
% [x,y,z,q] = READ_Q(filename)
% [r,q] = READ_Q(filename)
%
% Reads q=(\Delta c^{1/2})/c^{1/2} computed by the C++ codes (where c is
% a conductivity). For all but radsym and fast_radsym, it is computed at
% the nodes of a gmsh mesh.
%
% No q is computed for piecewise constant conductivities (pcc, fast_pcc,
% pcc_cmplx, fast_pcc_cmplx).
%
% 1) For real radially symmetric conductivities computed by radsym and
%    fast_radsym, the file is composed of a header where the first line is
%    the word "radial" and the second line is the number of radii in [0,1]
%    where q is computed. The file is then composed of two columns, the
%    first column gives the values of the the radii in [0,1] and the second
%    column gives the values q at the corresponding radii.
%
% 2) For other conductivities, the file is composed of a header which is
%    either a positive integer np only or a positive integer np followed on
%    the same line by the word "complex".
%    a) The first case corresponds to real q computed by the C++ codes
%       moments, fast_moments, eit, fast_eit.
%       After the header, the file is then composed of four columns. The
%       first three columns give cartesian coordinates of the points where
%       q is computed and the last column gives the values of q at the
%       corresponding points.
%    b) The second case corresponds to complex q computed by the C++ codes
%       moments_cmplx, fast_moments_cmplx, eit_cmplx, fast_eit_cmplx.
%       After the header, the file is then composed of five columns. The
%       first three columns give cartesian coordinates of the points where
%       q is computed and the last two columns give the values of the real
%       and imaginary part of q at the corresponding points.
%
% Rem: the function is exactly the same as read_conductivity.

% Open file
fid=fopen(filename,'r');

% Check error
if fid==-1
    error(['Error when opening file ',filename]);
end

% Read the first line and remove leading and trailing white space
firstline=strtrim(fgetl(fid));

if strcmp(firstline,'radial') % real radially symmetric conductivity
    np=fscanf(fid,'%u',1);fgetl(fid);
    rq=fscanf(fid,'%e',[2,np]).';
    varargout=num2cell(rq,1);
else % non-radial or complex conductivity
    [np xyzq]=strtok(firstline);
    np=str2double(np);
    xyzq=strtrim(xyzq);
    switch xyzq
        case '' % real conductivity
            xyzq=fscanf(fid,'%e',[4,np]).';
            varargout=num2cell(xyzq,1);
        case 'complex' % complex conductivity
            xyzq=fscanf(fid,'%e',[5,np]).';
            varargout=num2cell([xyzq(:,1:3),complex(xyzq(:,4),xyzq(:,5))],1);
        otherwise
            error(['Error while reading header of ',filename]);
    end
end

% Close file
fclose(fid);