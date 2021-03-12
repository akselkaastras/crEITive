#include "eit.h"

/*
  Interpolate at given nodes a function given on a uniform grid.
*/

void trilinear_interpolation(VecDouble &valnodes,MatDouble &nodes,Array3Double &val,const double &xmin,const double &xmax,const double &ymin,const double &ymax,const double &zmin,const double &zmax,const unsigned &nx,const unsigned &ny,const unsigned &nz)
{
  std::cout << "Trilinear interpolation..." << std::endl;

  unsigned n=nodes.size1();

  double dx=(xmax-xmin)/((double)nx-1.0);
  double dy=(ymax-ymin)/((double)ny-1.0);
  double dz=(zmax-zmin)/((double)nz-1.0);

  int ix,iy,iz;
  double tx,ty,tz;

  for (unsigned p=0;p<n;p++)
    {
      tx=(nodes(p,0)-xmin)/dx;
      ty=(nodes(p,1)-ymin)/dy;
      tz=(nodes(p,2)-zmin)/dz;
      //Cell
      ix=floor(tx);
      iy=floor(ty);
      iz=floor(tz);
      //Barycentric coordinates
      tx-=(double)ix;
      ty-=(double)iy;
      tz-=(double)iz;
      //Take sqrt(eps) here instead of eps
      //1e3*eps should be sufficient but in case...
      //If very close to a "left" face
      if (ix==-1&&1.0-tx<=sqrt(eps))
	{
	  ix=0;
	  tx=0;
	}
      if (iy==-1&&1.0-ty<=sqrt(eps))
	{
	  iy=0;
	  ty=0;
	}
      if (iz==-1&&1.0-tz<=sqrt(eps))
	{
	  iz=0;
	  tz=0;
	}
      //If very close to a "right" face
      if (ix==(int)nx-1&&tx<=sqrt(eps))
	{
	  ix=nx-2;
	  tx=1.0;
	}
      if (iy==(int)ny-1&&ty<=sqrt(eps))
	{
	  iy=ny-2;
	  ty=1.0;
	}
      if (iz==(int)nz-1&&tz<=sqrt(eps))
	{
	  iz=nz-2;
	  tz=1.0;
	}
      //Value
      valnodes(p)=(1.0-tx)*(1.0-ty)*(1.0-tz)*val[ix][iy][iz]+tx*(1.0-ty)*(1.0-tz)*val[ix+1][iy][iz]+(1.0-tx)*ty*(1.0-tz)*val[ix][iy+1][iz]+(1.0-tx)*(1.0-ty)*tz*val[ix][iy][iz+1]+tx*(1.0-ty)*tz*val[ix+1][iy][iz+1]+(1.0-tx)*ty*tz*val[ix][iy+1][iz+1]+tx*ty*(1.0-tz)*val[ix+1][iy+1][iz]+tx*ty*tz*val[ix+1][iy+1][iz+1];
    }
}
