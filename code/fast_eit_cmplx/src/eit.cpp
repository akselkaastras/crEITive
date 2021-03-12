#include "eit.h"
#include "read_commands.h"
#include "read_dnmap.h"
#include "gauss_legendre.h"
#include "phi_values.h"
#include "x_sphere.h"
#include "read_dnmap0.h"
#include "associated_legendre_functions.h"
#include "cos_m_phi.h"
#include "dirichlet_to_neumann_map_0.h"
#include "save_dnmap0.h"
#include "xi_grid.h"
#include "zeta_xi.h"
#include "qhat_zeta_exp.h"
#include "read_s0.h"
#include "slo_matrix.h"
#include "save_s0.h"
#include "qhat_zeta_0.h"
#include "qhat_zeta.h"
#include "save_psizeta.h"
#include "save_qhat.h"
#include "save_q.h"
#include "read_gmsh_mesh.h"
#include "trilinear_interpolation.h"
#include "ifft_calderon.h"
#include "whittaker_shannon_calderon.h"
#include "ifft_qhat.h"
#include "whittaker_shannon.h"
#include "schroedinger_matrix.h"
#include "schroedinger_rhs.h"
#include "solve_schroedinger.h"
#include "save_conductivity.h"

int main(int argc,char *argv[])
{
  //For CPU time
  double start = omp_get_wtime();

  //Commands filename
  string commandsfilename;

  //Log filename
  string logfilename;

  //Dirichlet-to-Neumann map filename
  string dnmapfilename;

  //Name of method used for the reconstruction
  string methodname;

  //Name of method for Inverse Fourier Transform
  string iftmethodname;

  //Name of method for the choice of zeta
  string zetamethodname;

  //Related to the quadrature points
  unsigned n;

  //Related to Fourier Transform grid
  unsigned ngrid;
  double pkappa;
  double truncrad;

  //Mesh filename to solve the Schrödinger equation
  //or to interpolate the reconstructed conductivity
  //if the Calderón approximation is chosen
  string meshfilename;

  //Name of the file where to save the CGO solutions traces
  string cgostracesfilename;
  //Name of the file where to save the reconstructed q
  string qfilename;
  //Name of the file where to save the reconstructed Fourier Transform of q
  string qhatfilename;

  //Name of the file where to save the reconstructed conductivity
  string conductivityfilename;

  //Optional
  //Dirichlet-to-Neumann map for a conductivity equal to 1
  string iodnmap0="";
  string dnmap0filename="";
  //Single Layer Operator Matrix
  string ios0="";
  string s0filename="";

  if (argc==1)
    {
      std::cout << "Enter commands filename: " << std::endl;
      std::cin >> commandsfilename;
      std::cout << "Enter log filename (none for display): " << std::endl;
      std::cin >> logfilename;
    }
  else if (argc==2)
    {
      commandsfilename=argv[1];
    }
  else if (argc==3)
    {
      commandsfilename=argv[1];
      logfilename=argv[2];
    }
  else
    {
      std::cerr << "Please give only 0, 1 or 2 arguments" << std::endl;
      exit(1);
    }

  if (logfilename=="none")
    {
      logfilename="";
    }

  std::ofstream trace_log (logfilename.c_str());
  std::streambuf *origcout,*origcerr;
  origcout=std::cout.rdbuf();
  origcerr=std::cerr.rdbuf();

  //If log file asked
  if (trace_log)
    {
      //Output stream
      //Connect stream buffers
      std::cout.rdbuf(trace_log.rdbuf());
      std::cerr.rdbuf(trace_log.rdbuf());
    }

  //Read commands
  read_commands(commandsfilename,ngrid,truncrad,iftmethodname,zetamethodname,pkappa,dnmapfilename,meshfilename,methodname,cgostracesfilename,qfilename,qhatfilename,conductivityfilename,iodnmap0,dnmap0filename,ios0,s0filename);
  
  //
  double scale = (double)(ngrid-1)*(double)(ngrid-1)*Pi/(2.0*(double)ngrid*(double)truncrad);
  std::cout << "the scale is " << scale << std::endl;
  
  //Dirichlet-to-Neumann map
  MatDoubleColMaj rdnmap;
  MatDoubleColMaj idnmap;
  read_dnmap(dnmapfilename,n,rdnmap,idnmap);

  //Angular quadrature points
  VecDouble tn (n+1);
  VecDouble wtn (n+1);
  gauss_legendre(tn,wtn);
  VecDouble phin (2*(n+1));
  phi_values(phin);

  //Weights
  wtn*=Pi/((double)(n+1));
  //Repeat them to give the weight for each quadrature point
  VecDouble wtnrep = rep_vec_1(wtn,2*(n+1));

  //Quadrature points
  MatDouble xn (2*square(n+1),3);
  x_sphere(xn,tn,phin);

  //Dirichlet-to-Neumann map for the conductivity 1
  MatDoubleColMaj dnmap0 (2*square(n+1),2*square(n+1));
  if (iodnmap0=="read") //Read it from file
    {
      read_dnmap0(dnmap0filename,dnmap0);
    }
  else //Compute it
    {
      //Associated Legendre functions
      MatDouble lfn (n+1,((n+1)*(n+2))/2);
      associated_legendre_functions(lfn,tn);
      //Cosines
      MatDouble cpn (2*(n+1),n+1);
      cos_m_phi(cpn,phin);
      dirichlet_to_neumann_map_0(dnmap0,lfn,cpn,wtn);
      if (iodnmap0=="write") //Save it
	{
	  save_dnmap0(n,dnmap0,dnmap0filename);
	}
    }

  //Difference of the dnmap
  rdnmap-=dnmap0;
  MatComplexColMaj dnmap = complexfunmat(rdnmap,idnmap);

  //kappa
  if (pkappa<1.0)
    {
      std::cout << "pkappa should be larger than 1" << std::endl;
      exit(1);
    }
  if (pkappa==1.0)
    {
      pkappa+=eps;
    }

  //xi, zeta
  MatDouble xi;
  xi_grid(xi,ngrid,iftmethodname,truncrad);
  unsigned ngrid3=pow(ngrid,3);
  VecDouble kappa (ngrid3);
  MatDouble k (ngrid3,3);
  MatDouble kT (ngrid3,3);
  zeta_xi(pkappa,xi,kappa,k,kT,zetamethodname);

  //CGO solutions traces
  MatComplexColMaj matpsizeta (xn.size1(),kappa.size());

  //Scattering Transform
  Array3Complex qhat (boost::extents[ngrid][ngrid][ngrid]);

  //Chosen reconstruction method

  //calderon or texp
  if (methodname=="calderon"||methodname=="texp")
    {
      qhat_zeta_exp(dnmap,xn,wtnrep,xi,kappa,k,kT,matpsizeta,qhat);
    }
  else if (methodname=="t0"||methodname=="t")
    {
      //Single Layer Operator Matrix
      MatDoubleColMaj S0 (2*(n+1)*(n+1),2*(n+1)*(n+1));
      if (ios0=="read") //Read it from file
	{
	  read_s0(s0filename,S0);
	}
      else //Compute it
	{
	  slo_matrix(S0,xn,wtnrep);
	  if (ios0=="write") //Save it
	    {
	      save_s0(n,S0,s0filename);
	    }
	}
      if (methodname=="t0")
	{
	  qhat_zeta_0(dnmap,rdnmap,idnmap,S0,xn,wtnrep,xi,kappa,k,kT,matpsizeta,qhat);
	}
      else if (methodname=="t")
	{
	  qhat_zeta(dnmap,rdnmap,idnmap,S0,xn,wtnrep,xi,kappa,k,kT,matpsizeta,qhat);
	}
    }

  //Save CGO solutions traces
  save_psizeta(matpsizeta,xi,kappa,k,kT,tn,phin,cgostracesfilename);
  
  //Save scattering transform
  save_qhat(qhat,xi,qhatfilename);

  //Read mesh
  unsigned nnodes, nquadnodes, nsurftriangles, nsurfnodes, ndof, ntetrahedra;
  VecUnsigned quadnodes, nodetoquad, surfnodes, nodetosurf, surftritet, surftriopp, dof, doftonode;
  MatUnsigned surftriangles, tetrahedra;
  VecDouble surftriarea, tetrahedravol;
  MatDouble nodes, surftribary, tetrahedrabary;
  VecMatDouble grad;
  read_gmsh_mesh(meshfilename,nnodes,nquadnodes,nsurftriangles,nsurfnodes,ndof,ntetrahedra,quadnodes,nodetoquad,surfnodes,nodetosurf,surftritet,surftriopp,dof,doftonode,surftriangles,tetrahedra,surftriarea,tetrahedravol,nodes,surftribary,tetrahedrabary,grad);

  VecComplex sig;

  if (methodname=="calderon")
    {
      sig.resize(nnodes);
      if (iftmethodname=="ifft") // IFFT method
	{
	  //Calderón formula with IFFT
	  ifft_calderon(xi,qhat);
	  //Trilinear interpolation
	  trilinear_interpolation(sig,nodes,qhat,-scale,scale,-scale,scale,-scale,scale,ngrid,ngrid,ngrid);
	}
      else if (iftmethodname=="ws") // Whittaker-Shannon method
	{
	  //Calderón formula with Whittaker-Shannon interpolation
	  whittaker_shannon_calderon(sig,nodes,xi,qhat);
	}
    }
  else if (methodname=="texp"||methodname=="t0"||methodname=="t")
    {
      sig.resize(ndof);
      VecComplex qbary (ntetrahedra);
      if (iftmethodname=="ifft") // IFFT method
	{
	  //IFFT
	  ifft_qhat(xi,qhat);
	  //Trilinear interpolation
	  trilinear_interpolation(sig,nodes,qhat,-scale,scale,-scale,scale,-scale,scale,ngrid,ngrid,ngrid);
	}
      else if (iftmethodname=="ws") // Whittaker-Shannon method
	{
	  //Whittaker-Shannon interpolation
	  whittaker_shannon(qbary,tetrahedrabary,xi,qhat);
	}
      // Save qbary
      save_q(qbary,tetrahedrabary,qfilename);

      //Schrödinger Matrix
      CompMatComplexColMaj MS (ndof,ndof,0);
      schroedinger_matrix(MS,qbary,nnodes,ntetrahedra,tetrahedra,tetrahedravol,dof,grad);

      //Schrödinger Right Hand Side
      sig.resize(ndof);
      schroedinger_rhs(sig,qbary,nnodes,ntetrahedra,tetrahedra,tetrahedravol,dof,grad);

      //Solve Schrödinger Equation
      solve_schroedinger(MS,sig,nnodes,ndof,nsurfnodes,doftonode,surfnodes);

      //Square to have the conductivity
      AutoVecOp(sig,square<dcomplex>);
    }

  //Save conductivity
  save_conductivity(nodes,sig,conductivityfilename);

  double end = omp_get_wtime();
  std::cout << "Elapsed time: " << end-start << "s" << std::endl;

  //If a trace log is opened, close it
  if (trace_log)
    {
      //Reconnect buffers to their original
      std::cout.rdbuf(origcout);
      std::cerr.rdbuf(origcerr);
      //Close file
      trace_log.close();
    }

  return 0;
}
