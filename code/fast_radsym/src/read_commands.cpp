#include "radsym.h"

/*
  Reads the commands file.
*/

void read_commands(string &commandsfilename,unsigned &nd,double &stepsize,MatDouble &sigma,string &dnmapfilename,string &dnmapeigenvaluesfilename,double &conductivitystepsize,string &conductivityfilename,double &qstepsize,string &qfilename)
{
  std::cout << "Reading the commands file..." << std::endl;

  std::ifstream commands_file(commandsfilename.c_str(),std::ios::in);

  if (!commands_file)
    {
      std::cerr << "Can't open file " << commandsfilename << std::endl;
      exit(1);
    }

  std::cout << "File " << commandsfilename << " successfully opened" << std::endl;

  string commands_line;
  string card;
  conductivityfilename="";
  qfilename="";
  qfilename="";
  unsigned nsigma=0;
  unsigned k=0;
  bool inconductivitycard=false;

  while(getline(commands_file,commands_line))
    {
      string card="";
      while (!commands_line.empty()) //Remove blank characters
	{
	  if (commands_line.substr(0,1)!=" ") //Write them in card
	    {
	      card=card+commands_line.substr(0,1);
	    }
	  commands_line.erase(0,1);
	}

      if (card=="#conductivity")
	{
	  k=nsigma;
	  nsigma+=1;
	  sigma.resize(nsigma,3);
	  inconductivitycard=true;
	}
      else if (card=="##radius"&&inconductivitycard)
	{
	  commands_file >> sigma(k,0);
	  getline(commands_file,commands_line);
	}
      else if (card=="##amplitude"&&inconductivitycard)
	{
	  commands_file >> sigma(k,1);
	  getline(commands_file,commands_line);
	}
      else if (card=="##innerextent"&&inconductivitycard)
	{
	  commands_file >> sigma(k,2);
	  getline(commands_file,commands_line);
	}
      else if (card=="#dnmap")
	{
	  inconductivitycard=false;
	  commands_file >> dnmapfilename;
	  getline(commands_file,commands_line);
	}
      else if (card=="#nd")
	{
	  inconductivitycard=false;
	  commands_file >> nd;
	  getline(commands_file,commands_line);
	}
      else if (card=="#stepsize")
	{
	  inconductivitycard=false;
	  commands_file >> stepsize;
	  getline(commands_file,commands_line);
	}
      else if (card=="#dnmapeigenvalues")
	{
	  inconductivitycard=false;
	  commands_file >> dnmapeigenvaluesfilename;
	  getline(commands_file,commands_line);
	}
      else if (card=="#conductivityonradius")
	{
	  inconductivitycard=false;
	  commands_file >> conductivitystepsize;
	  getline(commands_file,commands_line);
	  commands_file >> conductivityfilename;
	  getline(commands_file,commands_line);
	}
      else if (card=="#qonradius")
	{
	  inconductivitycard=false;
	  commands_file >> qstepsize;
	  getline(commands_file,commands_line);
	  commands_file >> qfilename;
	  getline(commands_file,commands_line);
	}
    }

  std::cout << "Parameters for the conductivity function" << std::endl;
  std::cout << "Number of components: " << nsigma << std::endl;
  for (k=0;k<nsigma;k++)
    {
      std::cout << "Component " << k+1 << std::endl;
      std::cout << "Radius: " << sigma(k,0) << std::endl;
      std::cout << "Amplitude: " << sigma(k,1) << std::endl;
      std::cout << "Inner extent: " << sigma(k,0) << std::endl;
    }
  std::cout << "Step size of the discretization to compute" << std::endl;
  std::cout << "the eigenvalues of the the Dirichlet-to-Neumann map" << std::endl;
  std::cout << "h=" << stepsize << std::endl;
  if (!dnmapeigenvaluesfilename.empty())
    {
      std::cout << "The eigenvalues of the Dirichlet-to-Neumann map" << std::endl;
      std::cout << "will be saved in the file:" << std::endl;
      std::cout << dnmapeigenvaluesfilename << " (in the directory dnmapeigenvalues/)" << std::endl;
    }
  if (!conductivityfilename.empty())
    {
      std::cout << "The conductivity function on a radius will be saved in the file:" << std::endl;
      std::cout << conductivityfilename << " (in the directory conductivity/)" << std::endl;
    }
  if (!qfilename.empty())
    {
      std::cout << "q function on a radius will be saved in the file:" << std::endl;
      std::cout << qfilename << " (in the directory q/)" << std::endl;
    }
  std::cout << "Parameters for Dirichlet-to-Neumann map" << std::endl;
  std::cout << "nd=" << nd << std::endl;
  std::cout << "The Dirichlet-to-Neumann map will be saved in the file:" << std::endl;
  std::cout << dnmapfilename << " (in the directory dnmap/)" << std::endl;
  if (!dnmapeigenvaluesfilename.empty())
    {
      std::cout << "The eigenvalues of the Dirichlet-to-Neumann map will be saved in the file:" << std::endl;
      std::cout << dnmapeigenvaluesfilename << " (in the directory dnmapeigenvalues/)" << std::endl;
    }

  commands_file.close();
}
