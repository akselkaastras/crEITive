#include "radsym.h"
#include "read_commands.h"
#include "conductivity_on_radius.h"
#include "q_on_radius.h"
#include "gauss_legendre.h"
#include "phi_values.h"
#include "cos_m_phi.h"
#include "associated_legendre_functions.h"
#include "dnmap_eigenvalues_pcc.h"
#include "dnmap_eigenvalues.h"
#include "dirichlet_to_neumann_map_radsym.h"
#include "save_conductivity.h"
#include "save_q.h"
#include "save_dnmap_eigenvalues.h"
#include "save_dnmap_radsym.h"

int main(int argc,char *argv[])
{
  //For execution time
  double start = omp_get_wtime();

  //Commands filename
  string commandsfilename;

  //Log filename
  string logfilename;

  //Name of the file where to save the Dirichlet-to-Neumann map
  string dnmapfilename;

  //Save Dirichlet-to-Neumann map eigenvalues
  string dnmapeigenvaluesfilename="";

  //Conductivity
  MatDouble sigma;

  //Save conductivity and/or q
  double conductivitystepsize;
  string conductivityfilename="";
  double qstepsize;
  string qfilename="";

  //Discretization
  double stepsize;

  //Related to Dirichlet-to-Neumann map
  unsigned nd;

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
  read_commands(commandsfilename,nd,stepsize,sigma,dnmapfilename,dnmapeigenvaluesfilename,conductivitystepsize,conductivityfilename,qstepsize,qfilename);

  //Compute and save conductivity on radius
  if (!conductivityfilename.empty())
    {
      unsigned nc=ceil(1.0/conductivitystepsize)+1;
      VecDouble r (nc);
      VecDouble valsigma (nc);
      //Compute
      conductivity_on_radius(r,valsigma,sigma);
      //Save
      save_conductivity(r,valsigma,conductivityfilename);
    }

  //Compute ans save q on radius
  if (!qfilename.empty())
    {
      unsigned nq=ceil(1.0/qstepsize)+1;
      VecDouble r (nq);
      VecDouble valq (nq);
      //Compute
      q_on_radius(r,valq,sigma);
      //Save
      save_q(r,valq,qfilename);
    }

  //Compute eigenvalues
  VecDouble eigenvalues (nd+1);
  VecDouble eigenvalueserror (nd+1);
  dnmap_eigenvalues(eigenvalues,eigenvalueserror,sigma,stepsize);

  //Save eigenvalues
  if (!dnmapeigenvaluesfilename.empty())
    {
      save_dnmap_eigenvalues(stepsize,eigenvalues,eigenvalueserror,dnmapeigenvaluesfilename);
    }

  //Compute Dirichlet-to-Neumann map
  MatDouble dnmap (2*(nd+1)*(nd+1),2*(nd+1)*(nd+1));

  //Angular quadrature points
  VecDouble tnd (nd+1);
  VecDouble wtnd (nd+1);
  gauss_legendre(tnd,wtnd);
  VecDouble phind (2*(nd+1));
  phi_values(phind);

  //Associated Legendre functions
  MatDouble lfnd (nd+1,((nd+1)*(nd+2))/2);
  associated_legendre_functions(lfnd,tnd,nd);

  //Cosines
  MatDouble cpnd (2*(nd+1),nd+1);
  cos_m_phi(cpnd,phind);
  dirichlet_to_neumann_map_radsym(dnmap,lfnd,cpnd,wtnd,eigenvalues);

  //Save Dirichlet-to-Neumann map
  save_dnmap_radsym(nd,dnmap,dnmapfilename);

  double end = omp_get_wtime();
  std::cout << "Elapsed time: " << end-start << "s" << std::endl;

  //If a trace log is opened, close it
  if (trace_log)
    {
      std::cout.rdbuf(origcout);
      std::cerr.rdbuf(origcerr);
      trace_log.close();
    }

  return 0;
}
