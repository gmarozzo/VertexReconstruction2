#include <Riostream.h>
#include "Retta.h"
#include <TMath.h>
#include <TRandom3.h>

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


void Retta::Rotate(Double_t th, Double_t ph, Double_t thp, Double_t php, Double_t *cd){
  double_t mr[3][3];

  mr[0][0] = -TMath::Sin(ph);
  mr[1][0] = TMath::Cos(ph);
  mr[2][0] = 0.;
  mr[0][1] = -TMath::Cos(ph)*TMath::Cos(th);
  mr[1][1] = -TMath::Cos(th)*TMath::Sin(ph);
  mr[2][1] = TMath::Sin(th);
  mr[0][2] = TMath::Sin(th)*TMath::Cos(ph);
  mr[1][2] = TMath::Sin(th)*TMath::Sin(ph);
  mr[2][2] = TMath::Cos(th);


  Double_t cdp[3];
  cdp[0] = TMath::Sin(thp)*TMath::Cos(php);
  cdp[1] = TMath::Sin(thp)*TMath::Sin(php);
  cdp[2] = TMath::Cos(thp);
  for(Int_t i = 0; i < 3; i++){
    cd[i]=0.;
    for(Int_t j = 0; j < 3; j++){
        cd[i] += mr[i][j]*cdp[j];
    }  
  }
}

void Retta::multiplescattering(){
  //gRandom = new TRandom3();
  double cd[3];
  double angs0=gRandom->Gaus(0,0.001);
  double angs1=(gRandom -> Rndm())*2.*M_PI;
  this->Rotate(this->GetTheta(),this->GetPhi(), angs0, angs1, cd);
  double thetanew = TMath::ACos(cd[2]);
  double phinew = atan2(cd[1], cd[0]);
  this->SetTheta(thetanew);
  this->SetPhi(phinew);
}