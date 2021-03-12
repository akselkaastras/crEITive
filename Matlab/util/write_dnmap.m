function write_dnmap(filename,DN)
% WRITE_DNMAP(filename,DNmap) writes a Dirichlet-to-Neumann possibly
% computed with Matlab or a Dirichlet-to-Neumann computed the C++ codes (
% moments, fast_moments, moments_cmplx, fast_moments_cmplx, pcc, fast_pcc,
% pcc_cmplx, fast_pcc_cmplx, radsym, fast_radsym) and modified with Matlab,
% typically to add noise, by first reading the Dirichlet-to-Neumann map
% with the function read_dnmap and then modifying it.
%
% The Dirichlet-to-Neumann map file is composed of a header which is either
% a positive integer nd only or a positive integer nd followed on the same
% line by the word "complex".
%
% 1) The first case corresponds to a real Dirichlet-to-Neumann map
%    matrix with 2*(nd+1)^2 quadrature points in the unit sphere. For the
%    After the header, the file is then composed of 2*(nd+1)^2 lines and
%    2*(nd+1)^2 columns giving the values of the Dirichlet-to-Neumann map
%    matrix: the j-th element of the (i+1)-th line is DN(i,j).
%
% 2) The second case corresponds to a complex Dirichlet-to-Neumann map
%    matrix with 2*(nd+1)^2 quadrature points in the unit sphere.
%    After the header, the file is then composed of 2*(nd+1)^2 lines and
%    2*(2*(nd+1)^2) columns giving the values of the Dirichlet-to-Neumann
%    map matrix: the (2*j)-th element of the (i+1)-th line is real(DN(i,j))
%    and the (2*j+1)-th element of the (i+1)-th line is imag(DN(i,j)).

% Open file
fid=fopen(filename,'w');

% Check error
if fid==-1
    error(['Error when opening file ',filename]);
end

% Write nd
np=size(DN,1);
if np~=size(DN,2);
    error('The matrix of the Dirichlet-to-Neumann map has to be square');
end

nd=sqrt(np/2)-1;
fprintf(fid,'%u',nd);

% If the matrix has complex at least a complex number, write "complex" on
% the next line and split each column of the matrix in two, the first for
% the real part, the second for the imaginary part.

DNc=whos('DN');
DNc=DNc.complex;

if DNc
    fprintf(fid,'%s\n',' complex');
    fclose(fid);
    DN=reshape([real(DN);imag(DN)],np,2*np);
    dlmwrite(filename,DN,'-append','precision','%.17e','delimiter',' ');
else
    fprintf(fid,'\n');
    fclose(fid);
    dlmwrite(filename,DN,'-append','precision','%.17e','delimiter',' ');
end
