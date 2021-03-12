#include "cgos.h"
#include "read_commands.h"
#include "read_dnmap.h"
#include "gauss_legendre.h"
#include "phi_values.h"
#include "x_sphere.h"
#include "read_dnmap0.h"
#include "cos_m_phi.h"
#include "associated_legendre_functions.h"
#include "dirichlet_to_neumann_map_0.h"
#include "save_dnmap0.h"
#include "read_s0.h"
#include "slo_matrix.h"
#include "save_s0.h"
#include "zeta_grid.h"
#include "cgos_traces_0.h"
#include "cgos_traces.h"
#include "save_psizeta.h"

int main(int argc,char *argv[])
{
  //For CPU time
  int thisprocess=RUSAGE_SELF;
  struct rusage usage;

  //Commands filename
  string commandsfilename;

  //Log filename
  string logfilename;

  //Dirichlet-to-Neumann map filename
  string dnmapfilename;

  //Name of method used to compute the CGO solutions traces
  string methodname;

  //Norm of zeta
  double zetanorm;

  //zetanorm/sqrt(2)*(kT+i*k), where k,kT \in \mathbb{R}^3
  // |k|=|kT|=1 and k.kT=0.
  //Number of discretization points of [0,2*pi) for the rotation around k
  //(number of discretization points for kT when k is chosen).
  unsigned nk;
  //Number of discretization points in [0,pi] for the elevation of k
  unsigned nt;
  //Number of discretization points in [0,2*pi) for the azimuth of k
  unsigned np;

  //Name of the file where to save the CGO solutions traces
  string cgostracesfilename;

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
  read_commands(commandsfilename,dnmapfilename,methodname,zetanorm,nk,nt,np,cgostracesfilename,iodnmap0,dnmap0filename,ios0,s0filename);

  //Read Dirichlet-to-Neumann map
  unsigned nd;
  MatDoubleColMaj dnmap;
  read_dnmap(dnmapfilename,nd,dnmap);

  //Angular quadrature points
  VecDouble tnd (nd+1);
  VecDouble wtnd (nd+1);
  gauss_legendre(tnd,wtnd);
  VecDouble phind (2*(nd+1));
  phi_values(phind);

  //Weights
  wtnd*=Pi/((double)(nd+1));
  //Repeat them to give the weight for each quadrature point
  VecDouble wtndrep = rep_vec_1(wtnd,2*(nd+1));

  //Quadrature points
  MatDouble xnd (2*square(nd+1),3);
  x_sphere(xnd,tnd,phind);

  //Dirichlet-to-Neumann map for the conductivity 1
  MatDoubleColMaj dnmap0 (2*square(nd+1),2*square(nd+1));
  if (iodnmap0=="read") //Read it from file
    {
      read_dnmap0(dnmap0filename,dnmap0);
    }
  else //Compute it
    {
      //Cosines
      MatDouble cpnd (2*(nd+1),nd+1);
      cos_m_phi(cpnd,phind);
      //Associated Legendre functions
      MatDouble lfnd (nd+1,((nd+1)*(nd+2))/2);
      associated_legendre_functions(lfnd,tnd);
      dirichlet_to_neumann_map_0(dnmap0,lfnd,cpnd,wtnd);
      if (iodnmap0=="write") //Save it
	{
	  save_dnmap0(nd,dnmap0,dnmap0filename);
	}
    }
  //Single Layer Operator Matrix
  MatDoubleColMaj S0 (2*square(nd+1),2*square(nd+1));
  if (ios0=="read") //Read it from file
    {
      read_s0(s0filename,S0);
    }
  else //Compute it
    {
      slo_matrix(S0,xnd,wtndrep);
      if (ios0=="write") //Save it
	{
	  save_s0(nd,S0,s0filename);
	}
    }

  //Difference of the dnmap
  dnmap-=dnmap0;

  //zeta
  double kappa=zetanorm/sqrt(2.0);
  VecDouble akT (nk);
  VecDouble theta (nt);
  VecDouble phi (np);
  unsigned ng=nk*((nt-1)*np+1);
  MatDouble k (ng,3);
  MatDouble kT (ng,3);
  zeta_grid(k,kT,akT,theta,phi);

  //CGO solutions traces
  MatComplexColMaj matpsizeta (2*square(nd+1),ng);

  //Chosen method

  if (methodname=="t0")
    {
      cgos_traces_0(dnmap,S0,xnd,wtndrep,kappa,k,kT,matpsizeta);
    }
  else if (methodname=="t")
    {
      cgos_traces(dnmap,S0,xnd,wtndrep,kappa,k,kT,matpsizeta);
    }

  //Save CGO solutions traces
  save_psizeta(matpsizeta,kappa,akT,theta,phi,tnd,phind,cgostracesfilename);

  getrusage(thisprocess,&usage);

  std::cout << "Elapsed time: " << usage.ru_utime.tv_sec << "." << usage.ru_utime.tv_usec << "s" << std::endl;

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
