function m = read_gmsh_mesh(filename)
% m = READ_GMSH_MESH(filename) reads a gmsh mesh file stored in the file
% filename.
%
% Allocates a mesh structure gmsh_mesh using the gmsh ASCII mesh file
% mesh_file.
%
% The gmsh ASCIII mesh file has the following form (remember that 15 is the
% element type for a point, 2 for a triangle, 4 for a 4-node tetrahedron):
%
% ...
% ...
% ...
% $Nodes
% 1 x1 y1 z1
% 2 x2 y2 z2
% ...
% $EndNodes
% $Elements
% 1 15 3 1 * * num_node
% 2 15 3 1 * * num_node
% ...
% ...
% p    2 3 2 * * triangle_nodes
% p+1  2 3 2 * * triangle_nodes
% ...
% ...
% q    4 3 3 * * tetrahedron_nodes
% q+1  4 3 3 * * tetrahedron_nodes
% $EndElements
%
% The forth integer is for the "physical group"
% 1) Special points, typically the quadrature points wich will be used for
%    the inverse code: allow to have directly nodes at these positions and
%    to list these particular nodes for a future use (legacy).
% 2) Triangles on the surface of the ball: points of these triangles are
%    not degrees of freedom
% 3) The ball
%
% m is a stucture with the following fields:
%
% * nnodes: an integer, the number of nodes
% * nodes: a (3 x nnodes) real matrix, each column for the coordinates of
%   the nodes
% * nquadnodes: an integer, the number of special nodes (group 1) 
% * quadnodes: a row integer vector (nquadnodes elements), the number of
%   each special node
% * quadtonode: a row integer vector (nnodes elements), for p,
%   quadtonode(p) is 0 if the node p is not a special node and it is its
%   number in the quadnodes list otherwise
% * nsurftriangles: an integer, the number of triangles on the surface of
%   the ball
% * surftriangles: a (3 x nsurftriangles) integer matrix, each column for
%   the nodes numbers of the surface triangles vertices
% * areasurftri: a real row vector (nsurftriangles elements), the area of
%   the surface triangles
% * nnodessurf: an integer, the number of nodes on the ball surface
% * nodessurf: a row integer vector (nnodesurf elements), the number of
%   each surface node
% * nodetosurf: a row integer vector (nnodes elements), for p,
%   nodetosurf(p) is 0 if the node p is not on the ball surface and it is
%   its number in the nodessurf list otherwise
% * surftritet: a row integer vector (nsurftriangles elements),
%   surftritet(p) is the number of the tetrahedron which contains the
%   triangle p
% * surftriopp: a row integer vector (nsurftriangles elements),
%   surftriopp(p) is the number of the opposite node to the triangle p in
%   the tetrahedron surftritet(p)
% * barysurftri: a (3 x nsurftriangles) real matrix, each column for
%   the barycenter coordinates of the surface triangles
% * ntetrahedra: an integer, the number of tetrahedra
% * tetrahedra: a (4 x ntetrahedra) integer matrix, each column for the
%   nodes numbers of the tetrahedra vertices
% * ndof: an integer, the number of degrees of freedom
% * dof: an integer row vector (nnodes elements), degree of freedom number
%   of each node (0 if a node is not a degree of freedom)
% * doftonode: an integer row vector (ndof elements), nodes numbers of the
%   different degrees of freedom 
% * barytetrahedra: a (3 x ntetrahedra) real matrix, each column for the
%   barycenter coordinates of the tetrahedra
% * voltetrahedra: a real row vector (ntetrahedra elements), the tetrahedra
%   volumes
% * grad: a (3 x 4 x ntetrahedra) real array, each page for each
%   tetrahedron, each column for the gradient of the hat function
%   associated to each vertex of the tetrahedron

% Open file
fid=fopen(filename,'r');

% Read number of nodes
nnodes='';
while ~strcmp(nnodes,'$Nodes')
    nnodes=fgetl(fid);
end
m.nnodes=fscanf(fid,'%d',1);

% Read nodes
m.nodes=fscanf(fid,'%d %g %g %g',[4,m.nnodes]);
m.nodes(1,:)=[];

% Read number of elements
m.ntetrahedra='';
while ~strcmp(m.ntetrahedra,'$Elements')
    m.ntetrahedra=fgetl(fid);
end
m.ntetrahedra=fscanf(fid,'%d',1);

% Read nodes in physical group 1 (quadrature points)
group=1;
m.quadnodes=[];
m.nodetoquad=zeros(1,m.nnodes);
while group==1
    position=ftell(fid);
    fidline=fscanf(fid,'%d %d %d %d %d %d %d',7);
    group=fidline(4);
    if group==1
        m.quadnodes = [m.quadnodes,fidline(end)];
        m.nodetoquad(fidline(end))=size(m.quadnodes,2);
    end
end
m.nquadnodes=size(m.quadnodes,2);
m.ntetrahedra=m.ntetrahedra-m.nquadnodes;

% Read triangles in physical group 2 (ball surface)
fseek(fid,position,'bof');
group=2;
m.surftriangles=[];
while group==2
    position=ftell(fid);
    fidline=fscanf(fid,'%d %d %d %d %d %d %d %d %d',9);
    group=fidline(4);
    if group==2
        m.surftriangles = [m.surftriangles,fidline(end-2:end)];
    end
end
m.nsurftriangles=size(m.surftriangles,2);
m.ntetrahedra=m.ntetrahedra-m.nsurftriangles;

% Read tetrahedra in physical group 3 (ball)
fseek(fid,position,'bof');
m.tetrahedra=fscanf(fid,'%d %d %d %d %d %d %d %d %d %d',[10,m.ntetrahedra]);
m.tetrahedra(1:6,:)=[];

% Little check
lastline=fscanf(fid,'%s');
if ~strcmp(lastline,'$EndElements')
    disp(' ');
    disp('End of file differs from what expected.');
    warning('MATLAB:read_gmsh_mesh','Check gmsh file.');
end

% Close file
fclose(fid);

% Post processing

% Compute degress of freedom
m.dof=ones(1,m.nnodes);
m.dof(reshape(m.surftriangles,1,3*m.nsurftriangles))=0;
m.ndof=nnz(m.dof);
m.doftonode=find(m.dof);
m.dof(logical(m.dof))=(1:m.ndof);

% Nodes on the ball surface
m.nnodessurf=m.nnodes-m.ndof;
m.nodessurf=find(m.dof(1:m.nnodes)==0);
m.nodetosurf=zeros(1,m.nnodes);
m.nodetosurf(m.nodessurf)=(1:m.nnodessurf);

% Vectorized computation of surftritet and surftriopp
logitet=m.dof(m.tetrahedra);
logitet=(logitet==0);
logitet3=(sum(logitet,1)==3);
logitet4=(sum(logitet,1)==4);
tet3=(1:m.ntetrahedra);
tet3=tet3(logitet3);
tet4=(1:m.ntetrahedra);
tet4=tet4(logitet4);
n3=nnz(logitet3);
n4=nnz(logitet4);
indtet3o4=[tet3,repmat(tet4,1,4)];

clear tet3 tet4;

sf3=repmat(logitet3,4,1);
sf3(~logitet)=false;
sf3=m.tetrahedra(sf3);
sf3=reshape(sf3,3,n3);

sf4=repmat(logitet4,4,1);
sf4=m.tetrahedra(sf4);
sf4=reshape(sf4,4,n4);
sf4=[sf4([1 2 3],:),sf4([2 3 4],:),sf4([3 4 1],:),sf4([4 1 2],:)];

clear logitet logitet3 logitet4 n3 n4;

sf=[sf3,sf4];
sf=sort(sf);

clear sf3 sf4;

tf=sort(m.surftriangles);

ff=[tf,sf];

clear tf sf;

[ff pp qq]=unique(ff.','rows');

clear ff;

m.surftritet=zeros(1,m.nsurftriangles);
m.surftritet=indtet3o4(pp(qq(1:m.nsurftriangles))-m.nsurftriangles);

cellsetxor=@(M,N)(setxor(M,N));
m.surftriopp=cellfun(cellsetxor,...
    num2cell(m.surftriangles,1),...
    num2cell(m.tetrahedra(:,m.surftritet),1)...
    );

% % Not vectorized computation of surftritet and surftriopp
% indsurftet=m.dof(m.tetrahedra);
% indsurftet=(indsurftet==0);
% indsurftet=any(indsurftet,1);
% indsurftet=find(indsurftet);
% surftet=m.tetrahedra(:,indsurftet);
% nsurftet=size(surftet,2);
% 
% m.surftritet=zeros(1,m.nsurftriangles);
% whichtet=false(3,nsurftet);
% clear nsurftet;
% for q=1:m.nsurftriangles
%     whichtet(1,:)=any(surftet==m.surftriangles(1,q),1);
%     whichtet(2,:)=any(surftet==m.surftriangles(2,q),1);
%     whichtet(3,:)=any(surftet==m.surftriangles(3,q),1);
%     m.surftritet(q)=indsurftet(all(whichtet,1));
%     m.surftriopp(q)=setxor(m.tetrahedra(:,m.surftritet(q)),m.surftriangles(:,q));
% end
% clear indsurftet surftet whichtet;

% Compute barycenters and areas of the surface triangles
m.barysurftri=reshape(m.surftriangles,1,3*m.nsurftriangles);
m.barysurftri=m.nodes(:,m.barysurftri);
m.barysurftri=reshape(m.barysurftri,[3,3,m.nsurftriangles]);
m.areasurftri=m.barysurftri;

% Barycenters: for a triangle "a,b,c" : ( a + b + c )/3
m.barysurftri=squeeze(sum(m.barysurftri,2))/3;

% Areas: for a triangle "a,b,c" :
% A = | ( b - a ) x ( c - a ) |/2
m.areasurftri=m.areasurftri(:,2:3,:)-...
    repmat(m.areasurftri(:,1,:),[1,2,1]);
m.areasurftri=sqrt(...
    sum(cross(m.areasurftri(:,1,:),m.areasurftri(:,2,:),1).^2,1)...
    )/2;

% Reshape m.areasurftri
m.areasurftri=reshape(m.areasurftri,[1,m.nsurftriangles]);

% Compute barycenters, volumes of tetrahedra, normal to faces
m.barytetrahedra=reshape(m.tetrahedra,1,4*m.ntetrahedra);
m.barytetrahedra=m.nodes(:,m.barytetrahedra);
m.barytetrahedra=reshape(m.barytetrahedra,[3,4,m.ntetrahedra]);
m.voltetrahedra=m.barytetrahedra;
m.grad=m.barytetrahedra;

% Barycenters: for a tetrahedron "a,b,c,d" : ( a + b + c + d )/4
m.barytetrahedra=squeeze(sum(m.barytetrahedra,2))/4;

% Volumes: for a tetrahedron "a,b,c,d" :
% V = | ( b - a ) . ( ( c - a ) x (d - a) ) |/6
m.voltetrahedra=m.voltetrahedra(:,2:4,:)-...
    repmat(m.voltetrahedra(:,1,:),[1,3,1]);
m.voltetrahedra=abs(dot(m.voltetrahedra(:,1,:),...
    cross(m.voltetrahedra(:,2,:),m.voltetrahedra(:,3,:),1),1))/6;

% Gradient of the hat functions in each tetrahedron
% The gradient of the hat function associated to a node "a" in a
% tetrahedron "a,b,c,d" is given by the inward normal to the
% opposite triangle "b,c,d" times the triangle surface, divived by
% three times the tetrahedron volume
% Hence, up to an eventual minus sign :
% \grad \hat{phi}_a = ( ( c - b ) x ( d - b ) )/( 6V )

m.grad=cross(...
    m.grad(:,[2 4 4 2],:)-m.grad(:,[3 3 1 1],:),...
    m.grad(:,[4 1 2 3],:)-m.grad(:,[2 3 4 1],:),...
    1)./(6*repmat(m.voltetrahedra,[3,4,1]));

% Reshape m.voltetrahedra
m.voltetrahedra=reshape(m.voltetrahedra,[1,m.ntetrahedra]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It seems that gmsh correctly orientes the tetrahedra, hence the last
% part is useless. If there is any problem :
% 1) Add check_inw=m.barytetrahedra; after m.grad=m.barytetrahedra;
% 2) Comment
%    m.voltetrahedra=reshape(m.voltetrahedra,[1,m.ntetrahedra]);
% 3) Uncomment the last part
% 4) Add at the end
%    m.voltetrahedra=reshape(m.voltetrahedra,[1,m.ntetrahedra]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inward normals to get the gradient, checking if it is the case
% % Perhaps not useful if the mesh is well oriented, check this point in
% % gmsh documentation!
% check_inw=cat(2,...
%     sum(...
%     dot(...
%     check_inw(:,[2 3 4],:)-repmat(check_inw(:,1,:),[1 3 1]),...
%     repmat(m.grad(:,1,:),[1 3 1])...
%     ,1)...
%     ,2),...
%     sum(...
%     dot(...
%     check_inw(:,[3 4 1],:)-repmat(check_inw(:,2,:),[1 3 1]),...
%     repmat(m.grad(:,2,:),[1 3 1])...
%     ,1)...
%     ,2),...
%     sum(...
%     dot(...
%     check_inw(:,[4 1 2],:)-repmat(check_inw(:,3,:),[1 3 1]),...
%     repmat(m.grad(:,3,:),[1 3 1])...
%     ,1)...
%     ,2),...
%     sum(...
%     dot(...
%     check_inw(:,[1 2 3],:)-repmat(check_inw(:,4,:),[1 3 1]),...
%     repmat(m.grad(:,4,:),[1 3 1])...
%     ,1)...
%     ,2)...
%     );
% check_inw=(check_inw>0);
% % nnz(check_inw)
% check_inw=repmat(check_inw,[3 1 1]);
% 
% m.grad(check_inw)=-m.grad(check_inw);