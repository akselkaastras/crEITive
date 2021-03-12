#include "eit.h"

/*
  Reads a gmsh mesh file, computes all necessary mesh values.

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

  * filename: gmsh filename
  * nnodes: an integer, the number of nodes
  * nodes: a (nnodes x 3) real matrix, each column for the coordinates of the nodes
  * nquadnodes: an integer, the number of special nodes (group 1) 
  * quadnodes: a row integer vector (nquadnodes elements), the number of each special node
  * nodetoquad: a row integer vector (nnodes elements), for p, quadtonode(p) is 0 if the node p is not a special node and it is its number in the quadnodes list otherwise
  * nsurftriangles: an integer, the number of triangles on the surface of the ball
  * surftriangles: a (nsurftriangles x 3) integer matrix, each column for the nodes numbers of the surface triangles vertices
  * surftriarea: a real row vector (nsurftriangles elements), the area of the surface triangles
  * nsurfnodes: an integer, the number of nodes on the ball surface
  * surfnodes: a row integer vector (nsurfnodes elements), the number of each surface node
  * nodetosurf: a row integer vector (nnodes elements), for p, nodetosurf(p) is nnodes if the node p is not on the ball surface and it is its number in the surfnodes list otherwise
  * surftritet: a row integer vector (nsurftriangles elements), surftritet(p) is the number of the tetrahedron which contains the surface triangle p
  * surftriopp: a row integer vector (nsurftriangles elements), surftriopp(p) is the number of the opposite node to the surface triangle p in the tetrahedron surftritet(p)
  * surftribary: a (nsurftriangles x 3) real matrix, each column for the barycenter coordinates of the surface triangles
  * ntetrahedra: an integer, the number of tetrahedra
  * tetrahedra: a (ntetrahedra x 4) integer matrix, each column for the nodes numbers of the tetrahedra vertices
  * ndof: an integer, the number of degrees of freedom
  * dof: an integer row vector (nnodes elements), degree of freedom number of each node (nnodes if a node is not a degree of freedom)
  * doftonode: an integer row vector (ndof elements), nodes numbers of the different degrees of freedom 
  * tetrahedrabary: a (3 x ntetrahedra) real matrix, each column for the barycenter coordinates of the tetrahedra
  * tetrahedravol: a real row vector (ntetrahedra elements), the tetrahedra volumes
  * grad: a  ntetrahedra vector of (4 x 3)  real matrices, each matrix for each tetrahedron, each column for the gradient of the hat function associated to each vertex of the tetrahedron
  */

void read_gmsh_mesh(string &filename,unsigned &nnodes,unsigned &nquadnodes,unsigned &nsurftriangles,unsigned &nsurfnodes,unsigned &ndof,unsigned &ntetrahedra,VecUnsigned &quadnodes,VecUnsigned &nodetoquad,VecUnsigned &surfnodes,VecUnsigned &nodetosurf,VecUnsigned &surftritet,VecUnsigned &surftriopp,VecUnsigned &dof,VecUnsigned &doftonode,MatUnsigned &surftriangles,MatUnsigned &tetrahedra,VecDouble &surftriarea,VecDouble &tetrahedravol,MatDouble &nodes,MatDouble &surftribary,MatDouble &tetrahedrabary,VecMatDouble &grad)
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


  std::cout << "Reading the gmsh file..." << std::endl;

  string gmsh_line;
  unsigned dummy;
  unsigned group;
  unsigned pos;

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

  // Reading number of elements
  gmsh_line="";
  while (gmsh_line!="$Elements")
    {
      getline(gmsh_file,gmsh_line);
    }
  gmsh_file >> ntetrahedra;
  getline(gmsh_file,gmsh_line);

  // Number of nodes in physical group 1 (quadrature points)
  group=1;
  nquadnodes=0;
  pos=gmsh_file.tellg();
  while (group==1)
    {
      gmsh_file >> dummy >> dummy >> dummy >> group;
      getline(gmsh_file,gmsh_line);
      if (group==1)
	{
	  nquadnodes++;
	}
    }

  // Reading nodes in physical group 1 (quadrature points)
  quadnodes.resize(nquadnodes);
  nodetoquad.resize(nnodes);
  for (unsigned i=0;i<nnodes;i++)
    {
      nodetoquad(i)=-1;
    }
  gmsh_file.seekg(pos,std::ios::beg);
  for (unsigned i=0;i<nquadnodes;i++)
    {
      gmsh_file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> quadnodes(i);
      getline(gmsh_file,gmsh_line);
      quadnodes(i)--;
      nodetoquad(quadnodes(i))=i;
    }

  ntetrahedra=ntetrahedra-nquadnodes;

  // Number of triangles in physical group 2 (ball surface)
  group=2;
  nsurftriangles=0;
  pos=gmsh_file.tellg();
  while (group==2)
    {
      gmsh_file >> dummy >> dummy >> dummy >> group;
      getline(gmsh_file,gmsh_line);
      if (group==2)
	{
	  nsurftriangles++;
	}
    }

  // Reading triangles in physical group 2 (ball surface)
  surftriangles.resize(nsurftriangles,3);
  gmsh_file.seekg(pos,std::ios::beg);
  for (unsigned i=0;i<nsurftriangles;i++)
    {
      gmsh_file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> surftriangles(i,0) >> surftriangles(i,1) >> surftriangles(i,2);
      getline(gmsh_file,gmsh_line);
      surftriangles(i,0)--;
      surftriangles(i,1)--;
      surftriangles(i,2)--;
    }

  ntetrahedra=ntetrahedra-nsurftriangles;

  // Reading tetrahedra in physical group 3 (ball)
  tetrahedra.resize(ntetrahedra,4);
  for (unsigned i=0;i<ntetrahedra;i++)
    {
      gmsh_file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> tetrahedra(i,0) >> tetrahedra(i,1) >> tetrahedra(i,2) >> tetrahedra(i,3);
      getline(gmsh_file,gmsh_line);
      tetrahedra(i,0)--;
      tetrahedra(i,1)--;
      tetrahedra(i,2)--;
      tetrahedra(i,3)--;
    }

  // Little check
  getline(gmsh_file,gmsh_line);
  if (gmsh_line!="$EndElements")
    {
      std::cerr << "End of file differs from what expected, check gmsh file." << std::endl;
      exit(1);
    }

  // Closing file
  gmsh_file.close();

  // Post processing

  // Degress of freedom and nodes on the ball surface
  dof.resize(nnodes);
  doftonode.resize(nnodes);
  nodetosurf.resize(nnodes);
  surfnodes.resize(nnodes);
  for (unsigned i=0;i<nnodes;i++)
    {
      dof(i)=i;
    }
  for (unsigned i=0;i<nsurftriangles;i++)
    {
      dof(surftriangles(i,0))=nnodes;
      dof(surftriangles(i,1))=nnodes;
      dof(surftriangles(i,2))=nnodes;
    }
  nsurfnodes=0;
  for (unsigned i=0;i<nnodes;i++)
    {
      nodetosurf(i)=nnodes;
    }
  for (unsigned i=0;i<nnodes;i++)
    {
      if (dof(i)==nnodes)
	{
	  surfnodes(nsurfnodes)=i;
	  nodetosurf(nsurfnodes)=i;
	  nsurfnodes++;
	}
      else
	{
	  dof(i)-=nsurfnodes;
	  doftonode(dof(i))=i;
	}
    }
  surfnodes.resize(nsurfnodes,true);
  ndof=nnodes-nsurfnodes;
  doftonode.resize(ndof,true);

  // Computation of surftritet and surftriopp
  surftritet.resize(nsurftriangles);
  surftriopp.resize(nsurftriangles);
  for (unsigned i=0;i<ntetrahedra;i++)
    {
      unsigned node0=tetrahedra(i,0);
      unsigned node1=tetrahedra(i,1);
      unsigned node2=tetrahedra(i,2);
      unsigned node3=tetrahedra(i,3);
      bool bnode0=(dof(node0)==nnodes);
      bool bnode1=(dof(node1)==nnodes);
      bool bnode2=(dof(node2)==nnodes);
      bool bnode3=(dof(node3)==nnodes);

      if (bnode0&&bnode1&&bnode2)
	{
	  for (unsigned ii=0;ii<nsurftriangles;ii++)
	    {
	      unsigned snode0=surftriangles(ii,0);
	      unsigned snode1=surftriangles(ii,1);
	      unsigned snode2=surftriangles(ii,2);
	      if ((snode0==node0||snode0==node1||snode0==node2)&&(snode1==node0||snode1==node1||snode1==node2)&&(snode2==node0||snode2==node1||snode2==node2))
		{
		  surftritet(ii)=i;
		  surftriopp(ii)=node3;
		}
	    }
	}

      if (bnode0&&bnode1&&bnode3)
	{
	  for (unsigned ii=0;ii<nsurftriangles;ii++)
	    {
	      unsigned snode0=surftriangles(ii,0);
	      unsigned snode1=surftriangles(ii,1);
	      unsigned snode2=surftriangles(ii,2);
	      if ((snode0==node0||snode0==node1||snode0==node3)&&(snode1==node0||snode1==node1||snode1==node3)&&(snode2==node0||snode2==node1||snode2==node3))
		{
		  surftritet(ii)=i;
		  surftriopp(ii)=node2;
		}
	    }
	}

      if (bnode0&&bnode2&&bnode3)
	{
	  for (unsigned ii=0;ii<nsurftriangles;ii++)
	    {
	      unsigned snode0=surftriangles(ii,0);
	      unsigned snode1=surftriangles(ii,1);
	      unsigned snode2=surftriangles(ii,2);
	      if ((snode0==node0||snode0==node2||snode0==node3)&&(snode1==node0||snode1==node2||snode1==node3)&&(snode2==node0||snode2==node2||snode2==node3))
		{
		  surftritet(ii)=i;
		  surftriopp(ii)=node1;
		}
	    }
	}

      if (bnode1&&bnode2&&bnode3)
	{
	  for (unsigned ii=0;ii<nsurftriangles;ii++)
	    {
	      unsigned snode0=surftriangles(ii,0);
	      unsigned snode1=surftriangles(ii,1);
	      unsigned snode2=surftriangles(ii,2);
	      if ((snode0==node1||snode0==node2||snode0==node3)&&(snode1==node1||snode1==node2||snode1==node3)&&(snode2==node1||snode2==node2||snode2==node3))
		{
		  surftritet(ii)=i;
		  surftriopp(ii)=node0;
		}
	    }
	}
    }

  // Barycenters of the surface tringles
  // For a triangle "a,b,c": G = ( a + b + c )/3
  surftribary.resize(nsurftriangles,3);
  for (unsigned i=0;i<nsurftriangles;i++)
    {
      surftribary(i,0)=(nodes(surftriangles(i,0),0)+nodes(surftriangles(i,1),0)+nodes(surftriangles(i,2),0))/3.0;
      surftribary(i,1)=(nodes(surftriangles(i,0),1)+nodes(surftriangles(i,1),1)+nodes(surftriangles(i,2),1))/3.0;
      surftribary(i,2)=(nodes(surftriangles(i,0),2)+nodes(surftriangles(i,1),2)+nodes(surftriangles(i,2),2))/3.0;
    }

  // Areas of the surface triangles
  // For a triangle "a,b,c": A = | ( b - a ) x ( c - a ) |/2
  surftriarea.resize(nsurftriangles);
  for (unsigned i=0;i<nsurftriangles;i++)
    {
      VecDouble bma (3);
      bma(0)=nodes(surftriangles(i,1),0)-nodes(surftriangles(i,0),0);
      bma(1)=nodes(surftriangles(i,1),1)-nodes(surftriangles(i,0),1);
      bma(2)=nodes(surftriangles(i,1),2)-nodes(surftriangles(i,0),2);
      VecDouble cma (3);
      cma(0)=nodes(surftriangles(i,2),0)-nodes(surftriangles(i,0),0);
      cma(1)=nodes(surftriangles(i,2),1)-nodes(surftriangles(i,0),1);
      cma(2)=nodes(surftriangles(i,2),2)-nodes(surftriangles(i,0),2);
      surftriarea(i)=sqrt(square(bma(1)*cma(2)-bma(2)*cma(1))+square(bma(2)*cma(0)-bma(0)*cma(2))+square(bma(0)*cma(1)-bma(1)*cma(0)))/2.0;
    }

  // Barycenters of the tetrahedra
  // For a tetrahedron "a,b,c,d": G = ( a + b + c + d )/4
  tetrahedrabary.resize(ntetrahedra,3);
  for (unsigned i=0;i<ntetrahedra;i++)
    {
      tetrahedrabary(i,0)=(nodes(tetrahedra(i,0),0)+nodes(tetrahedra(i,1),0)+nodes(tetrahedra(i,2),0)+nodes(tetrahedra(i,3),0))/4.0;
      tetrahedrabary(i,1)=(nodes(tetrahedra(i,0),1)+nodes(tetrahedra(i,1),1)+nodes(tetrahedra(i,2),1)+nodes(tetrahedra(i,3),1))/4.0;
      tetrahedrabary(i,2)=(nodes(tetrahedra(i,0),2)+nodes(tetrahedra(i,1),2)+nodes(tetrahedra(i,2),2)+nodes(tetrahedra(i,3),2))/4.0;
    }

  // Volumes of the tetrahedra
  // For a tetrahedron "a,b,c,d": V = | ( b - a ) . ( ( c - a ) x (d - a) ) |/6
  tetrahedravol.resize(ntetrahedra);
  for (unsigned i=0;i<ntetrahedra;i++)
    {
      VecDouble bma (3);
      bma(0)=nodes(tetrahedra(i,1),0)-nodes(tetrahedra(i,0),0);
      bma(1)=nodes(tetrahedra(i,1),1)-nodes(tetrahedra(i,0),1);
      bma(2)=nodes(tetrahedra(i,1),2)-nodes(tetrahedra(i,0),2);
      VecDouble cma (3);
      cma(0)=nodes(tetrahedra(i,2),0)-nodes(tetrahedra(i,0),0);
      cma(1)=nodes(tetrahedra(i,2),1)-nodes(tetrahedra(i,0),1);
      cma(2)=nodes(tetrahedra(i,2),2)-nodes(tetrahedra(i,0),2);
      VecDouble dma (3);
      dma(0)=nodes(tetrahedra(i,3),0)-nodes(tetrahedra(i,0),0);
      dma(1)=nodes(tetrahedra(i,3),1)-nodes(tetrahedra(i,0),1);
      dma(2)=nodes(tetrahedra(i,3),2)-nodes(tetrahedra(i,0),2);

      tetrahedravol(i)=fabs(bma(0)*(cma(1)*dma(2)-cma(2)*dma(1))+bma(1)*(cma(2)*dma(0)-cma(0)*dma(2))+bma(2)*(cma(0)*dma(1)-cma(1)*dma(0)))/6.0;
    }

  // Gradient of the hat functions in each tetrahedron
  /*
    The gradient of the hat function associated to a node "a" in a
    tetrahedron "a,b,c,d" is given by the inward normal to the
    opposite triangle "b,c,d" times the triangle surface, divived by
    three times the tetrahedron volume.
    Hence, up to an possible minus sign:
    \grad \hat{phi}_a = ( ( c - b ) x ( d - b ) )/( 6V )
  */

  grad.resize(ntetrahedra);
  for (unsigned i=0;i<ntetrahedra;i++)
    {
      MatDouble cmb (4,3);
      cmb(0,0)=nodes(tetrahedra(i,1),0)-nodes(tetrahedra(i,2),0);
      cmb(0,1)=nodes(tetrahedra(i,1),1)-nodes(tetrahedra(i,2),1);
      cmb(0,2)=nodes(tetrahedra(i,1),2)-nodes(tetrahedra(i,2),2);
      cmb(1,0)=nodes(tetrahedra(i,3),0)-nodes(tetrahedra(i,2),0);
      cmb(1,1)=nodes(tetrahedra(i,3),1)-nodes(tetrahedra(i,2),1);
      cmb(1,2)=nodes(tetrahedra(i,3),2)-nodes(tetrahedra(i,2),2);
      cmb(2,0)=nodes(tetrahedra(i,3),0)-nodes(tetrahedra(i,0),0);
      cmb(2,1)=nodes(tetrahedra(i,3),1)-nodes(tetrahedra(i,0),1);
      cmb(2,2)=nodes(tetrahedra(i,3),2)-nodes(tetrahedra(i,0),2);
      cmb(3,0)=nodes(tetrahedra(i,1),0)-nodes(tetrahedra(i,0),0);
      cmb(3,1)=nodes(tetrahedra(i,1),1)-nodes(tetrahedra(i,0),1);
      cmb(3,2)=nodes(tetrahedra(i,1),2)-nodes(tetrahedra(i,0),2);
      MatDouble dmb (4,3);
      dmb(0,0)=nodes(tetrahedra(i,3),0)-nodes(tetrahedra(i,1),0);
      dmb(0,1)=nodes(tetrahedra(i,3),1)-nodes(tetrahedra(i,1),1);
      dmb(0,2)=nodes(tetrahedra(i,3),2)-nodes(tetrahedra(i,1),2);
      dmb(1,0)=nodes(tetrahedra(i,0),0)-nodes(tetrahedra(i,2),0);
      dmb(1,1)=nodes(tetrahedra(i,0),1)-nodes(tetrahedra(i,2),1);
      dmb(1,2)=nodes(tetrahedra(i,0),2)-nodes(tetrahedra(i,2),2);
      dmb(2,0)=nodes(tetrahedra(i,1),0)-nodes(tetrahedra(i,3),0);
      dmb(2,1)=nodes(tetrahedra(i,1),1)-nodes(tetrahedra(i,3),1);
      dmb(2,2)=nodes(tetrahedra(i,1),2)-nodes(tetrahedra(i,3),2);
      dmb(3,0)=nodes(tetrahedra(i,2),0)-nodes(tetrahedra(i,0),0);
      dmb(3,1)=nodes(tetrahedra(i,2),1)-nodes(tetrahedra(i,0),1);
      dmb(3,2)=nodes(tetrahedra(i,2),2)-nodes(tetrahedra(i,0),2);
      MatDouble crosscmbdmb (4,3);
      crosscmbdmb(0,0)=cmb(0,1)*dmb(0,2)-cmb(0,2)*dmb(0,1);
      crosscmbdmb(0,1)=cmb(0,2)*dmb(0,0)-cmb(0,0)*dmb(0,2);
      crosscmbdmb(0,2)=cmb(0,0)*dmb(0,1)-cmb(0,1)*dmb(0,0);
      crosscmbdmb(1,0)=cmb(1,1)*dmb(1,2)-cmb(1,2)*dmb(1,1);
      crosscmbdmb(1,1)=cmb(1,2)*dmb(1,0)-cmb(1,0)*dmb(1,2);
      crosscmbdmb(1,2)=cmb(1,0)*dmb(1,1)-cmb(1,1)*dmb(1,0);
      crosscmbdmb(2,0)=cmb(2,1)*dmb(2,2)-cmb(2,2)*dmb(2,1);
      crosscmbdmb(2,1)=cmb(2,2)*dmb(2,0)-cmb(2,0)*dmb(2,2);
      crosscmbdmb(2,2)=cmb(2,0)*dmb(2,1)-cmb(2,1)*dmb(2,0);
      crosscmbdmb(3,0)=cmb(3,1)*dmb(3,2)-cmb(3,2)*dmb(3,1);
      crosscmbdmb(3,1)=cmb(3,2)*dmb(3,0)-cmb(3,0)*dmb(3,2);
      crosscmbdmb(3,2)=cmb(3,0)*dmb(3,1)-cmb(3,1)*dmb(3,0);

      grad(i)=crosscmbdmb/(6.0*tetrahedravol(i));
    }

  /*
    It seems that gmsh correctly orientes the tetrahedra.
    Check it in the documentation or check the orientation
    in the code with some lines.
  */
}
