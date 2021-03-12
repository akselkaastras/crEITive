#ifndef TRILINEAR_INTERPOLATION
#define TRILINEAR_INTERPOLATION

void trilinear_interpolation(VecComplex &valnodes,MatDouble &nodes,Array3Complex &val,const double &xmin,const double &xmax,const double &ymin,const double &ymax,const double &zmin,const double &zmax,const unsigned &nx,const unsigned &ny,const unsigned &nz);

#endif // TRILINEAR_INTERPOLATION
