#include "eit.h"

/*
  Computes the grid points xi where the scattering transform will be computed. 
*/

void xi_grid(MatDouble &xi,unsigned &ngrid,string &iftmethodname, double &ximaxonmo)
{

  //For the IFFT method
  if (iftmethodname=="ifft")
    {
      unsigned ngrid2=square(ngrid);
      unsigned ngrid3=ngrid*ngrid2;
      xi.resize(ngrid3,3);
      int ng=(int)ngrid-1;
      ximaxonmo = (double) ximaxonmo/((double)ng);
      //ximaxonmo= Pi*(double)ng/(2.0*(double)ngrid);
      //std::cout << ximaxonmo << "  is truncrad divided by ngrid-1" << std::endl;
      //double ximaxonmo=Pi*(double)ng/(2.0*(double)ngrid);
      VecDouble xiloop (3);

      for (int ii=-ng;ii<=ng;ii+=2)
	{
	  xiloop(0)=(double)ii*ximaxonmo;
	  for (int jj=-ng;jj<=ng;jj+=2)
	    {
	      xiloop(1)=(double)jj*ximaxonmo;
	      for (int kk=-ng;kk<=ng;kk+=2)
		{
		  xiloop(2)=(double)kk*ximaxonmo;

		  boostublas::row(xi,ngrid2*(ii+ng)/2+ngrid*(jj+ng)/2+(kk+ng)/2)=xiloop;
		}
	    }
	}
  
    }
  //For the WS method
  else if (iftmethodname=="ws")
    {
      int ng=(int)ngrid;
      ngrid=2*ngrid+1;
      int twongpone=(int)ngrid;
      int twongpone2=square(twongpone);
      int twongpone3=twongpone*twongpone2;
      xi.resize(twongpone3,3);
      VecDouble xiloop (3);

      for (int ii=-ng;ii<=ng;ii++)
	{
	  xiloop(0)=Pi*(double)ii;
	  for (int jj=-ng;jj<=ng;jj++)
	    {
	      xiloop(1)=Pi*(double)jj;
	      for (int kk=-ng;kk<=ng;kk++)
		{
		  xiloop(2)=Pi*(double)kk;

		  boostublas::row(xi,twongpone2*(ii+ng)+twongpone*(jj+ng)+(kk+ng))=xiloop;
		}
	    }
	}
    }
}
