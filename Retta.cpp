#include <Riostream.h>
#include "Retta.h"
#include <TMath.h>

using namespace TMath; 

ClassImp(Retta)

//default constructor
Retta::Retta():TObject(),
  theta(0.),
  phi(0.),
  c1(0.),
  c2(0.),
  c3(0.),
  X0(0.),
  Y0(0.),
  Z0(0.){
  }

//standard constructor
Retta::Retta(double *X, double *ang):TObject(),
  theta(ang[0]),
  phi(ang[1]),
  c1(Sin(ang[0])*Cos(ang[1])),
  c2(Sin(ang[0])*Sin(ang[1])),
  c3(Cos(ang[0])),
  X0(X[0]),
  Y0(X[1]),
  Z0(X[2]){
  }

					      
double Retta :: parameter (double R){  //valutazione parametro t
  
  double prod = X0*c1+Y0*c2;
  double C = c1*c1+c2*c2;
  double Delta = (prod*prod)-C*((X0*X0)+(Y0*Y0)-(R*R));
  if (Delta<0){
    cout<<"Errore, discriminante negativo."<<endl;	
    return 0;}
  double t;
  double t1 = (-prod+sqrt(Delta))/C;
  double t2 = (-prod-sqrt(Delta))/C;
  if(t1>0) t=t1;
  else t=t2;
  return t;
}


void Retta :: intpoint(double *X, double R){  //valutazione punto di intersezione
	
	double t =  this->parameter(R);

	X[0]=X0+c1*t;
	X[1]=Y0+c2*t;
	X[2]=Z0+c3*t;
	
}
