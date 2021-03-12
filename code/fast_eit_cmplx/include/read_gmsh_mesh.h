#ifndef READ_GMSH_MESH
#define READ_GMSH_MESH

void read_gmsh_mesh(string &filename,unsigned &nnodes,unsigned &nquadnodes,unsigned &nsurftriangles,unsigned &nsurfnodes,unsigned &ndof,unsigned &ntetrahedra,VecUnsigned &quadnodes,VecUnsigned &nodetoquad,VecUnsigned &surfnodes,VecUnsigned &nodetosurf,VecUnsigned &surftritet,VecUnsigned &surftriopp,VecUnsigned &dof,VecUnsigned &doftonode,MatUnsigned &surftriangles,MatUnsigned &tetrahedra,VecDouble &surftriarea,VecDouble &tetrahedravol,MatDouble &nodes,MatDouble &surftribary,MatDouble &tetrahedrabary,VecMatDouble &grad);

#endif // READ_GMSH_MESH
