#include "eit.h"

/*
  Reads the commands file.
*/

void read_commands(string &commandsfilename,unsigned &ngrid, double &truncrad,string &iftmethodname,string &zetamethodname,double &pkappa,string &dnmapfilename,string &meshfilename,string &methodname,string &cgostracesfilename,string &qfilename,string &qhatfilename,string &conductivityfilename,string &iodnmap0,string &dnmap0filename,string &ios0,string &s0filename)
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
		  if (card=="ngrid")
		    {
		      commands_file >> ngrid;
		      getline(commands_file,commands_line);
		    }
		  if (card=="truncrad")
		    {
		      commands_file >> truncrad;
		      getline(commands_file,commands_line);
		    }			
		  if (card=="ift")
		    {
		      commands_file >> iftmethodname;
		      getline(commands_file,commands_line);
		    }
		  if (card=="zeta")
		    {
		      commands_file >> zetamethodname;
		      getline(commands_file,commands_line);
		    }
		  if (card=="pkappa")
		    {
		      commands_file >> pkappa;
		      getline(commands_file,commands_line);
		    }
		  if (card=="method")
		    {
		      commands_file >> methodname;
		      getline(commands_file,commands_line);
		    }
		  if (card=="mesh")
		    {
		      commands_file >> meshfilename;
		      getline(commands_file,commands_line);
		    }
		  if (card=="q")
		    {
		      commands_file >> qfilename;
		      getline(commands_file,commands_line);
		    }
		  if (card=="qhat")
		    {
		      commands_file >> qhatfilename;
		      getline(commands_file,commands_line);
		    }
		  if (card=="conductivity")
		    {
		      commands_file >> conductivityfilename;
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

  if (methodname=="calderon")
    {
      std::cout << "Reconstruction with Calder�n approximation" << std::endl;
    }
  if (methodname=="texp")
    {
      std::cout << "Reconstruction with texp approximation" << std::endl;
    }
  if (methodname=="t0")
    {
      std::cout << "Reconstruction with t0 approximation" << std::endl;
    }
  if (methodname=="t")
    {
      std::cout << "Reconstruction with Nachman algorithm" << std::endl;
    }
  std::cout << "Truncation with ideal low-pass filter" << std::endl;
  std::cout << "truncation radius =" << (double) truncrad << std::endl;
  if (iftmethodname=="ifft")
    {
      std::cout << "Inverse Fourier Transform will be computed by IFFT" << std::endl;
      std::cout << "Number of grid points for the scattering transform ngrid^3:" << std::endl;
      std::cout << "ngrid=" << ngrid << std::endl;
	  std::cout << "truncation radius=" << truncrad << std::endl;
    }
  if (iftmethodname=="ws")
    {
      std::cout << "Inverse Fourier Transform will be computed with Whittaker-Shannon formula" << std::endl;
      std::cout << "Number of grid points for the scattering transform (2*ngrid+1)^3:" << std::endl;
      std::cout << "ngrid=" << ngrid << std::endl;
    }

  if (zetamethodname=="fixed")
    {
      std::cout << "Fixed norm zeta:" << std::endl;
      if (iftmethodname=="ifft")
	{
	  std::cout << "|zeta|=" << pkappa*truncrad/(sqrt(2.0)) << std::endl; 
	}
      if (iftmethodname=="ws")
	{
	  std::cout << "|zeta|=" << pkappa*sqrt(3)*Pi*(double)ngrid/sqrt(2.0) << std::endl;
	}
    }
  if (zetamethodname=="proportional")
    {
      std::cout << "Proportional norm zeta:" << std::endl;
      std::cout << "sqrt(2)*|zeta|/|xi|=" << pkappa << std::endl;
    }

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

  if (methodname=="t0"||methodname=="t")
    {
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
    }

  std::cout << "Mesh filename for interpolation:" << std::endl;
  std::cout << meshfilename << std::endl;
  std::cout << "The trace of the CGO solutions" << std::endl;
  std::cout << "will be saved in the file:" << std::endl;
  std::cout << cgostracesfilename << " (in the directory cgostraces/)" << std::endl;
  std::cout << "The reconstructed Fourier Transform of q will be saved in the file:" << std::endl;
  std::cout << qhatfilename << " (in the directory qhat/)" << std::endl;
  std::cout << "The reconstructed conductivity will be saved in the file:" << std::endl;
  std::cout << conductivityfilename << " (in the directory conductivity/)" << std::endl;

  commands_file.close();
}
