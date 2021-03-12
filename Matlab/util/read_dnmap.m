function [DN,cmplx,nd] = read_dnmap(filename)
% DN = READ_DNMAP(filename) reads a Dirichlet-to-Neumann map computed by
% the C++ codes:
% moments, fast_moments, moments_cmplx, fast_moments_cmplx, pcc, fast_pcc,
% pcc_cmplx, fast_pcc_cmplx, radsym, fast_radsym.
%
% The Dirichlet-to-Neumann map file is composed of a header which is either
% a positive integer nd only or a positive integer nd followed on the same
% line by the word "complex".
%
% 1) The first case corresponds to a real Dirichlet-to-Neumann map
%    matrix with 2*(nd+1)^2 quadrature points in the unit sphere. This case
%    is for the C++ codes moments, fast_moments, pcc, fast_pcc,
%    radsym, fast_radsym.
%    After the header, the file is then composed of 2*(nd+1)^2 lines and
%    2*(nd+1)^2 columns giving the values of the Dirichlet-to-Neumann map
%    matrix: the j-th element of the (i+1)-th line is DN(i,j).
%
% 2) The second case corresponds to a complex Dirichlet-to-Neumann map
%    matrix with 2*(nd+1)^2 quadrature points in the unit sphere. This case
%    is for the C++ codes moments_cmplx, fast_moments_cmplx,
%    pcc_cmplx, fast_pcc_cmplx.
%    After the header, the file is then composed of 2*(nd+1)^2 lines and
%    2*(2*(nd+1)^2) columns giving the values of the Dirichlet-to-Neumann
%    map matrix: the (2*j)-th element of the (i+1)-th line is real(DN(i,j))
%    and the (2*j+1)-th element of the (i+1)-th line is imag(DN(i,j)).

% Open file
fid=fopen(filename,'r');

% Check error
if fid==-1
    error(['Error when opening file ',filename]);
end

% Read nd
nd=fscanf(fid,'%u',1);
np=2*(nd+1)^2;

% Read the end of the line and remove leading and trailing white space
rori=strtrim(fgetl(fid));

% Close file
fclose(fid);

switch rori
    case '' % real matrix
        DN=dlmread(filename,' ',1,0);
    case 'complex' % complex matrix
        DN=dlmread(filename,' ',1,0); 
        DN=complex(DN(:,1:2:end),DN(:,2:2:end));
    otherwise
        error(['Error while reading header of ',filename]);
end

if ~(np==size(DN,1)&&np==size(DN,2))
    error('The matrix of the Dirichlet-to-Neumann map has to be square');
end

cmplx = strcmp(rori,'complex');