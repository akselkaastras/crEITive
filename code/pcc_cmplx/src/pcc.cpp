#include "pcc.h"
#include "read_commands.h"
#include "read_gmsh_mesh_nodes.h"
#include "conductivity_at_nodes.h"
#include "save_conductivity.h"
#include "gauss_legendre.h"
#include "associated_legendre_functions.h"
#include "grad_solid_harmonics.h"
#include "phi_values.h"
#include "exp_i_m_phi.h"
#include "kp_matrix.h"
#include "kpcross_matrix.h"
#include "kbarp_matrix.h"
#include "kp_kbarp_matrix.h"
#include "rhs_matrix.h"
#include "rhs_block_matrix.h"
#include "solve_pcc.h"
#include "dev_to_val.h"
#include "block_dev_to_val.h"
#include "normal_derivatives.h"
#include "dirichlet_to_neumann_map.h"
#include "save_dnmap.h"

int main(int argc,char *argv[])
{
  //For execution time
  clock_t start = clock();

  //Commands filename
  string commandsfilename;

  //Log filename
  string logfilename;

  //Name of the file where to save the Dirichlet-to-Neumann map
  string dnmapfilename;

  //Conductivity
  VecConductivity sigma;

  //Save conductivity
  string conductivitymeshfilename="";
  string conductivityfilename="";
  int nnodes;
  MatDouble nodes;
  VecComplex sigmanodes;

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
  read_commands(commandsfilename,nd,sigma,dnmapfilename,conductivitymeshfilename,conductivityfilename);

  //Save conductivity
  if (!conductivitymeshfilename.empty()&&!conductivityfilename.empty())
    {
      read_gmsh_mesh_nodes(conductivitymeshfilename,nnodes,nodes);
      conductivity_at_nodes(sigmanodes,nodes,sigma);
      save_conductivity(nodes,sigmanodes,conductivityfilename);
    }

  //Matrix
  unsigned matsize=0;
  for (unsigned k=0;k<sigma.size();k++)
    {
      matsize+=pow(sigma(k).GetMaxDegree()+1,2);
    }
  MatComplex M (matsize,matsize);
  kp_kbarp_matrix(M,sigma);

  //Right hand side
  MatComplex rhs (matsize,(nd+1)*(nd+1));
  rhs_block_matrix(rhs,sigma,nd);

  //Solve linear system
  solve_pcc(M,rhs);

  //Values of densities
  MatComplex valdens (2*matsize,(nd+1)*(nd+1));
  block_dev_to_val(valdens,rhs,sigma);
  rhs.resize(0,0);

  //Quadrature points on the unit sphere
  VecDouble td (nd+1);
  VecDouble costd (nd+1);
  VecDouble sintd (nd+1);
  VecDouble pd (2*(nd+1));
  VecDouble cospd (2*(nd+1));
  VecDouble sinpd (2*(nd+1));
  VecDouble alphad (nd+1);
  MatDouble xd (2*(nd+1)*(nd+1),3);
  gauss_legendre(costd,alphad);
  alphad*=Pi/((double)(nd+1));
  std::transform(costd.begin(),costd.end(),td.begin(),[](double c) -> double { return acos(c); });
  std::transform(td.begin(),td.end(),sintd.begin(),[](double c) -> double { return sin(c); });
  phi_values(pd);
  std::transform(pd.begin(),pd.end(),cospd.begin(),[](double c) -> double { return cos(c); });
  std::transform(pd.begin(),pd.end(),sinpd.begin(),[](double c) -> double { return sin(c); });
  unsigned nid;
  for (unsigned ipd=0;ipd<2*(nd+1);ipd++)
    {
      for (unsigned itd=0;itd<nd+1;itd++)
	{
	  nid=itd+(nd+1)*ipd;
	  xd(nid,0)=sintd(itd)*cospd(ipd);xd(nid,1)=sintd(itd)*sinpd(ipd);xd(nid,2)=costd(itd);
	}
    }

  //Normal derivatives
  MatComplex ND (2*(nd+1)*(nd+1),(nd+1)*(nd+1));
  normal_derivatives(ND,valdens,sigma,xd);

  //Dirichlet-to-Neumann map
  MatComplex DNmap (2*(nd+1)*(nd+1),2*(nd+1)*(nd+1));
  MatDouble lfd (nd+1,((nd+1)*(nd+2))/2);
  MatComplex epd (2*(nd+1),2*nd+1);
  associated_legendre_functions(lfd,costd,nd);
  exp_i_m_phi(epd,pd);
  dirichlet_to_neumann_map(DNmap,ND,lfd,epd,alphad);

  //Save Dirichlet-to-Neumann map
  save_dnmap(nd,DNmap,dnmapfilename);

  std::cout << "Elapsed time: " << ((double)clock()-start)/CLOCKS_PER_SEC << "s" << std::endl;

  //If a trace log is opened, close it
  if (trace_log)
    {
      std::cout.rdbuf(origcout);
      std::cerr.rdbuf(origcerr);
      trace_log.close();
    }

  return 0;
}
