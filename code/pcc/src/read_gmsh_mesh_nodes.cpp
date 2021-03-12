#include "pcc.h"

/*
  Reads the nodes of a gmsh mesh file.

  The gmsh ASCII mesh file has the following form (remember that 15 is the element type for a point, 2 for a triangle, 4 for a 4-node tetrahedron):

  ...
  ...
  ...
  $Nodes
  1 x1 y1 z1
  2 x2 y2 z2
  ...
  $EndNodes
  $Elements
  1 15 3 1 * * num_node
  2 15 3 1 * * num_node
  ...
  ...
  p    2 3 2 * * triangle_nodes
  p+1  2 3 2 * * triangle_nodes
  ...
  ...
  q    4 3 3 * * tetrahedron_nodes
  q+1  4 3 3 * * tetrahedron_nodes
  $EndElements

  The forth integer is for the "physical group".
  1) Special points, possibly the quadrature points used for the inverse code: allow to have directly nodes at these positions and to list these particular nodes for a future use.
  2) Triangles on the surface of the ball (points of these triangles are not degrees of freedom).
  3) The ball.

  The following variables are returned by the function:

  * nnodes: an integer, the number of nodes
  * nodes: a (nnodes x 3) real matrix, each column for the coordinates of the nodes
  */

void read_gmsh_mesh_nodes(string &filename,int &nnodes,MatDouble &nodes)
{

  string fullfilename="mesh/"+filename;

  // Opening file
  std::ifstream gmsh_file(fullfilename.c_str(),std::ios::in);

  if (!gmsh_file)
    {
      std::cerr << "Can't open file " << fullfilename << std::endl;
      exit(1);
    }

  std::cout << "File " << fullfilename << " successfully opened" << std::endl;


  std::cout << "Reading the gmsh file nodes..." << std::endl;

  string gmsh_line;
  int dummy;

  // Reading number of nodes
  gmsh_line="";
  while (gmsh_line!="$Nodes")
    {
      getline(gmsh_file,gmsh_line);
    }
  gmsh_file >> nnodes;
  getline(gmsh_file,gmsh_line);

  // Reading nodes
  nodes.resize(nnodes,3);
  for (unsigned i=0;i<nnodes;i++)
    {
      gmsh_file >> dummy >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
      getline(gmsh_file,gmsh_line);
    }

  // Closing file
  gmsh_file.close();
}
