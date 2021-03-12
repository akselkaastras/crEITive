#include "cgos.h"

/*
  Computes the values of zeta where the CGO solutions traces are computed.
  zeta=kappa*(kT+i*k), with kappa=|zeta|/sqrt(2), |kT|=|k|=1, kT.k=0.
*/

void zeta_grid(MatDouble &k,MatDouble &kT,VecDouble &akT,VecDouble &theta,VecDouble &phi)
{
  std::cout << "Computing the values of zeta..." << std::endl;

  unsigned nk=akT.size();
  unsigned nt=theta.size();
  unsigned np=phi.size();

  double tpionk=2.0*Pi/(double)nk;
  double piont=Pi/(double)nt;
  double pionp=Pi/(double)np;

  unsigned ntm1=nt-1;
  VecDouble thetam (ntm1);

  std::copy(CountItDouble0,CountItDouble0+akT.size(),akT.begin());
  std::copy(CountItDouble1,CountItDouble1+thetam.size(),thetam.begin());
  std::copy(CountItDouble0,CountItDouble0+phi.size(),phi.begin());
  akT*=tpionk;
  thetam*=piont;
  phi*=pionp;

  theta(0)=0.0;
  boostublas::subrange(theta,1,nt)=thetam;

  VecDouble ckT=VecOp(akT,[](double c) -> double { return cos(c); });
  VecDouble skT=VecOp(akT,[](double c) -> double { return sin(c); });
  VecDouble ct=VecOp(thetam,[](double c) -> double { return cos(c); });
  VecDouble st=VecOp(thetam,[](double c) -> double { return sin(c); });
  VecDouble cp=VecOp(phi,[](double c) -> double { return cos(c); });
  VecDouble sp=VecOp(phi,[](double c) -> double { return sin(c); });

  unsigned ngw0=nk*ntm1*np;
  MatDouble kw0 (ngw0,3);
  MatDouble kTw0 (ngw0,3);

  //For theta!=0
  //kw0
  boostublas::column(kw0,0)=rep_vec_2(kron_prod_vec(cp,st),nk);
  boostublas::column(kw0,1)=rep_vec_2(kron_prod_vec(sp,st),nk);
  boostublas::column(kw0,2)=rep_vec_2(rep_vec_1(ct,np),nk);

  //kTw0
  boostublas::column(kTw0,0)=kron_prod_vec(kron_prod_vec(cp,ct),ckT)+rep_vec(kron_prod_vec((VecDouble)(-sp),skT),nk,ntm1);
  boostublas::column(kTw0,1)=kron_prod_vec((VecDouble)kron_prod_vec(sp,ct),ckT)+rep_vec(kron_prod_vec(cp,skT),nk,ntm1);
  boostublas::column(kTw0,2)=rep_vec_1(kron_prod_vec((VecDouble)(-st),ckT),np);

  //k
  boostublas::matrix_vector_slice<MatDouble>(k,boostublas::slice(0,1,nk),boostublas::slice(0,0,nk))=ZeroVecDouble(nk);
  boostublas::matrix_vector_slice<MatDouble>(k,boostublas::slice(nk,1,ngw0),boostublas::slice(0,0,ngw0))=boostublas::column(kw0,0);
  boostublas::matrix_vector_slice<MatDouble>(k,boostublas::slice(0,1,nk),boostublas::slice(1,0,nk))=ZeroVecDouble(nk);
  boostublas::matrix_vector_slice<MatDouble>(k,boostublas::slice(nk,1,ngw0),boostublas::slice(1,0,ngw0))=boostublas::column(kw0,1);
  boostublas::matrix_vector_slice<MatDouble>(k,boostublas::slice(0,1,nk),boostublas::slice(2,0,nk))=ScalVecDouble(nk,1.0);
  boostublas::matrix_vector_slice<MatDouble>(k,boostublas::slice(nk,1,ngw0),boostublas::slice(2,0,ngw0))=boostublas::column(kw0,2);

  //kT
  boostublas::matrix_vector_slice<MatDouble>(kT,boostublas::slice(0,1,nk),boostublas::slice(0,0,nk))=ckT;
  boostublas::matrix_vector_slice<MatDouble>(kT,boostublas::slice(nk,1,ngw0),boostublas::slice(0,0,ngw0))=boostublas::column(kTw0,0);
  boostublas::matrix_vector_slice<MatDouble>(kT,boostublas::slice(0,1,nk),boostublas::slice(1,0,nk))=skT;
  boostublas::matrix_vector_slice<MatDouble>(kT,boostublas::slice(nk,1,ngw0),boostublas::slice(1,0,ngw0))=boostublas::column(kTw0,1);
  boostublas::matrix_vector_slice<MatDouble>(kT,boostublas::slice(0,1,nk),boostublas::slice(2,0,nk))=ZeroVecDouble(nk);
  boostublas::matrix_vector_slice<MatDouble>(kT,boostublas::slice(nk,1,ngw0),boostublas::slice(2,0,ngw0))=boostublas::column(kTw0,2);
}
