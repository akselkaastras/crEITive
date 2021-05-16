#include "pcc.h"

/*
  Reads the commands file.
*/

void read_commands(string &commandsfilename,unsigned &nd,VecConductivity &sigma,string &dnmapfilename,string &conductivitymeshfilename,string &conductivityfilename)
{
  std::cout << "Reading the commands file and" << std::endl;
  std::cout << "and initializing surfaces data..." << std::endl;

  std::ifstream commands_file(commandsfilename.c_str(),std::ios::in);

  if (!commands_file)
    {
      std::cerr << "Can't open file " << commandsfilename << std::endl;
      exit(1);
    }

  std::cout << "File " << commandsfilename << " successfully opened" << std::endl;

  double x,y,z,r,theta,phi,ampli,normvec,width;
  unsigned n;
  string commands_line;
  string card;
  string coordsys;
  conductivitymeshfilename="";
  conductivityfilename="";
  unsigned nsigma=0;
  unsigned k=0;
  unsigned naxis;
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
	  sigma.resize(nsigma);
	  inconductivitycard=true;
	}
      else if (card=="##center"&&inconductivitycard)
	{
	  commands_file >> coordsys;
	  getline(commands_file,commands_line);
	  if (coordsys=="cartesian")
	    {
	      commands_file >> x;
	      commands_file >> y;
	      commands_file >> z;
	      getline(commands_file,commands_line);
	      sigma(k).SetCenterCoord(0,x);
	      sigma(k).SetCenterCoord(1,y);
	      sigma(k).SetCenterCoord(2,z);
	    }
	  else if (coordsys=="spherical")
	    {
	      commands_file >> r;
	      commands_file >> theta;
	      commands_file >> phi;
	      getline(commands_file,commands_line);
	      theta*=Pi/180.0;
	      phi*=Pi/180.0;
	      sigma(k).SetCenterCoord(0,r*sin(theta)*cos(phi));
	      sigma(k).SetCenterCoord(1,r*sin(theta)*sin(phi));
	      sigma(k).SetCenterCoord(2,r*cos(theta));
	    }
	}
      else if (card=="##radius"&&inconductivitycard)
	{
	  commands_file >> r;
	  getline(commands_file,commands_line);
	  sigma(k).SetAxes(IdMatDouble (3));
	  sigma(k).SetRadii(ScalVecDouble (3,r));
	}
      else if ((card=="##axis1"||card=="##axis2"||card=="##axis3")&&inconductivitycard)
	{
	  naxis=std::atoi(card.substr(6,1).c_str())-1;
	  commands_file >> coordsys;
	  getline(commands_file,commands_line);
	  if (coordsys=="cartesian")
	    {
	      commands_file >> x;
	      commands_file >> y;
	      commands_file >> z;
	      getline(commands_file,commands_line);
	      normvec=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	      sigma(k).SetAxisCoord(naxis,0,x/normvec);
	      sigma(k).SetAxisCoord(naxis,1,y/normvec);
	      sigma(k).SetAxisCoord(naxis,2,z/normvec);
	    }
	  else if (coordsys=="spherical")
	    {
	      commands_file >> theta;
	      commands_file >> phi;
	      getline(commands_file,commands_line);
	      theta*=Pi/180.0;
	      phi*=Pi/180.0;
	      sigma(k).SetAxisCoord(naxis,0,sin(theta)*cos(phi));
	      sigma(k).SetAxisCoord(naxis,1,sin(theta)*sin(phi));
	      sigma(k).SetAxisCoord(naxis,2,cos(theta));
	    }
	}
      else if ((card=="##radius1"||card=="##radius2"||card=="##radius3")&&inconductivitycard)
	{
	  naxis=std::atoi(card.substr(8,1).c_str())-1;
	  commands_file >> r;
	  getline(commands_file,commands_line);
	  sigma(k).SetRadius(naxis,r);
	}
	else if (card=="##width"&&inconductivitycard)
	{
	  commands_file >> width;
	  getline(commands_file,commands_line);
	  sigma(k).SetWidth(width);
	}
      else if (card=="##amplitude"&&inconductivitycard)
	{
	  commands_file >> ampli;
	  getline(commands_file,commands_line);
	  sigma(k).SetAmplitude(ampli);
	}
      else if (card=="##n"&&inconductivitycard)
	{
	  commands_file >> n;
	  getline(commands_file,commands_line);
	  sigma(k).SetMaxDegree(n);
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
      else if (card=="#conductivityonmesh")
	{
	  inconductivitycard=false;
	  commands_file >> conductivitymeshfilename;
	  getline(commands_file,commands_line);
	  commands_file >> conductivityfilename;
	  getline(commands_file,commands_line);
	}
    }

  std::cout << "Parameters for the conductivity function" << std::endl;
  std::cout << "Number of components: " << nsigma << std::endl;
  for (k=0;k<nsigma;k++)
    {
      std::cout << "Component " << k+1 << std::endl;
      sigma(k).Display();
    }
  if (!conductivitymeshfilename.empty()&&!conductivityfilename.empty())
    {
      std::cout << "The conductivity function at the points of" << std::endl;
      std::cout << "the mesh " <<  conductivitymeshfilename << " will be saved in the file:" << std::endl;
      std::cout << conductivityfilename << " (in the directory conductivity/)" << std::endl;
    }
  std::cout << "Parameters for Dirichlet-to-Neumann map" << std::endl;
  std::cout << "nd=" << nd << std::endl;
  std::cout << "The Dirichlet-to-Neumann map will be saved in the file:" << std::endl;
  std::cout << dnmapfilename << " (in the directory dnmap/)" << std::endl;

  commands_file.close();

  std::cout << "Initializing conductivity data..." << std::endl;
  for (unsigned k=0;k<sigma.size();k++)
    {
      std::cout << "Component " << k+1 << std::endl;
      sigma(k).InitializeData();
    }
}
