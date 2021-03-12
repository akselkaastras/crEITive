#include "cgos.h"

/*
  Reads the commands file.
*/

void read_commands(string &commandsfilename,string &dnmapfilename,string &methodname,double &zetanorm,unsigned &nk,unsigned &nt,unsigned &np,string &cgostracesfilename,string &iodnmap0,string &dnmap0filename,string &ios0,string &s0filename)
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

  //Optional cards initialized as empty
  iodnmap0="";
  dnmap0filename="";
  ios0="";
  s0filename="";

  while(getline(commands_file,commands_line))
    {

      if (!commands_line.empty()) //If line is not empty
	{
	  card="";
	  while (!commands_line.empty()) //Remove blank characters
	    {
	      if (commands_line.substr(0,1)!=" ") //Write them in card
		{
		  card=card+commands_line.substr(0,1);
		}
	      commands_line.erase(0,1);
	    }

	  if (card.substr(0,1)=="#") //If first character is #
	    {
	      card.erase(0,1);

	      if (!card.empty())
		{
		  if (card=="dnmap")
		    {
		      commands_file >> dnmapfilename;
		      getline(commands_file,commands_line);
		    }
		  if (card=="method")
		    {
		      commands_file >> methodname;
		      getline(commands_file,commands_line);
		    }
		  if (card=="zetanorm")
		    {
		      commands_file >> zetanorm;
		      getline(commands_file,commands_line);
		    }
		  if (card=="nk")
		    {
		      commands_file >> nk;
		      getline(commands_file,commands_line);
		    }
		  if (card=="nt")
		    {
		      commands_file >> nt;
		      getline(commands_file,commands_line);
		    }
		  if (card=="np")
		    {
		      commands_file >> np;
		      getline(commands_file,commands_line);
		    }
		  if (card=="cgostraces")
		    {
		      commands_file >> cgostracesfilename;
		      getline(commands_file,commands_line);
		    }
		  if (card=="dnmap0")
		    {
		      commands_file >> iodnmap0;
		      getline(commands_file,commands_line);
		      commands_file >> dnmap0filename;
		      getline(commands_file,commands_line);
		    }
		  if (card=="s0")
		    {
		      commands_file >> ios0;
		      getline(commands_file,commands_line);
		      commands_file >> s0filename;
		      getline(commands_file,commands_line);
		    }
		}
	    }
	}
    }

  std::cout << "Dirichlet-to-Neumann map filename:" << std::endl;
  std::cout << dnmapfilename << std::endl;

  if (methodname=="t0")
    {
      std::cout << "CGO solutions traces computed with t0 approximation" << std::endl;
    }
  if (methodname=="t")
    {
      std::cout << "CGO solutions traces computed with Nachman algorithm" << std::endl;
    }

  std::cout << "Zeta norm:" << std::endl;
  std::cout << "|zeta|=" << zetanorm << std::endl;

  std::cout << "zeta=|zeta|/sqrt(2)*(kT+i*k)" << zetanorm << std::endl;
  std::cout << "Number of discretization points in [0,2*pi) for the rotation around k" << std::endl;
  std::cout << "(number of discretization points for kT when k is fixed):" << std::endl;
  std::cout << "nk=" << nk << std::endl;
  std::cout << "Number of discretization points in [0,pi] for the elevation of k:" << std::endl;
  std::cout << "nt=" << nt << std::endl;
  std::cout << "Number of discretization points in [0,2*pi) for the azimuth of k:" << std::endl;
  std::cout << "np=" << np << std::endl;

  if (iodnmap0.empty())
    {
      std::cout << "The Dirichlet-to-Neumann map for a conductivity equal to 1" << std::endl;
      std::cout << "will be computed and not saved." << std::endl;
    }
  else
    {
      if (iodnmap0=="write")
	{
	  std::cout << "The Dirichlet-to-Neumann map for a conductivity equal to 1" << std::endl;
	  std::cout << "will be saved in the file:" << std::endl;
	  std::cout << dnmap0filename << " (in the directory dnmap0/)" << std::endl;
	}
      else if (iodnmap0=="read")
	{
	  std::cout << "The Dirichlet-to-Neumann map for a conductivity equal to 1" << std::endl;
	  std::cout << "will be read from the file:" << std::endl;
	  std::cout << dnmap0filename << " (in the directory dnmap0/)" << std::endl;
	}
    }

  if (ios0.empty())
    {
      std::cout << "The Single Layer Operator Matrix" << std::endl;
      std::cout << "will be computed and not saved." << std::endl;
    }
  else
    {
      if (ios0=="write")
	{
	  std::cout << "The Single Layer Operator Matrix" << std::endl;
	  std::cout << "will be saved in the file:" << std::endl;
	  std::cout << s0filename << " (in the directory s0/)" << std::endl;
	}
      else if (ios0=="read")
	{
	  std::cout << "The Single Layer Operator Matrix" << std::endl;
	  std::cout << "will be read from the file:" << std::endl;
	  std::cout << s0filename << " (in the directory s0/)" << std::endl;
	}
    }

  std::cout << "The traces of the CGO solutions will be saved in the file:" << std::endl;
  std::cout << cgostracesfilename << " (in the directory cgostraces/)" << std::endl;

  commands_file.close();
}
