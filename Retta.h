#ifndef RETTA_H
#define RETTA_H

#include "TObject.h"

class Retta: public TObject {


private:
  double theta, phi;
  double c1, c2, c3;
  double X0,Y0,Z0;
  void SetTheta(double th){theta=th;}
  void SetPhi(double ph){phi=ph;}
  void Rotate(Double_t th, Double_t ph, Double_t thp, Double_t php, Double_t *cd);
	
public:
  Retta();
  Retta(double *X, double *ang);
  double parameter(double R);
  void intpoint(double *X, double R);
  double GetTheta(){return theta;}
  double GetPhi(){return phi;}
  void multiplescattering();
};

#endif
