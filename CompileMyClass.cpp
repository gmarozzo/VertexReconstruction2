#include <TString.h>
#include <TSystem.h>
void compilemyclass(TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  gSystem->CompileMacro("MyRandom3.cpp",opt.Data());
  gSystem->CompileMacro("Vertex2.cpp",opt.Data());
  gSystem->CompileMacro("Hit2.cpp",opt.Data());
  gSystem->CompileMacro("Retta.cpp",opt.Data());
}
