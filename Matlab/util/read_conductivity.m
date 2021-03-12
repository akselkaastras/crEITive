function varargout = read_conductivity(filename)
% [x,y,z,c] = READ_CONDUCTIVITY(filename)
% [r,c] = READ_CONDUCTIVITY(filename)
%
% Reads a conductivity computed by the C++ codes. For all but radsym and
% fast_radsym, the conducvitity is computed at the nodes of a gmsh mesh.
%
% 1) For real radially symmetric conductivities computed by radsym and
%    fast_radsym, the file is composed of a header where the first line is
%    the word "radial" and the second line is the number of radii in [0,1]
%    where the radially symmetric conductivity is computed. The file is
%    then composed of two columns, the first column gives the values of the
%    the radii in [0,1] and the second column gives the values of the
%    conductivity at the corresponding radii.
%
% 2) For other conductivities, the file is composed of a header which is
%    either a positive integer np only or a positive integer np followed on
%    the same line by the word "complex".
%    a) The first case corresponds to real conductivities computed by the
%       C++ codes moments, fast_moments, pcc, fast_pcc, eit, fast_eit.
%       After the header, the file is then composed of four columns. The
%       first three columns give cartesian coordinates of the points where
%       the conductivity is computed and the last column gives the values
%       of the conductivity at the corresponding points.
%    b) The second case corresponds to complex conductivities computed by
%       the C++ codes moments_cmplx, fast_moments_cmplx, pcc_cmplx,
%       fast_pcc_cmplx, eit_cmplx, fast_eit_cmplx.
%       After the header, the file is then composed of five columns. The
%       first three columns give cartesian coordinates of the points where
%       the conductivity is computed and the last two columns give the
%       values of the real and imaginary part of the conductivity at the
%       corresponding points.
%
% Rem: the function is exactly the same as read_q.

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
    rc=fscanf(fid,'%e',[2,np]).';
    varargout=num2cell(rc,1);
else % non-radially-symmetric or complex conductivity
    [np xyzc]=strtok(firstline);
    np=str2double(np);
    xyzc=strtrim(xyzc);
    switch xyzc
        case ''
            xyzc=fscanf(fid,'%e',[4,np]).';
            varargout=num2cell(xyzc,1);
        case 'complex'
            xyzc=fscanf(fid,'%e',[5,np]).';
            varargout=num2cell([xyzc(:,1:3),complex(xyzc(:,4),xyzc(:,5))],1);
        otherwise
            error(['Error while reading header of ',filename]);
    end  
end

% Close file
fclose(fid);