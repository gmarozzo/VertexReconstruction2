#ifndef RETTA_H
#define RETTA_H

#include "TObject.h"

class Retta: public TObject {


private:
  double theta, phi;
  double c1, c2, c3;
  double X0,Y0,Z0;
	
public:
  Retta();
  Retta(double *X, double *ang);
  double parameter(double R);
  void intpoint(double *X, double R);
};

#endif
