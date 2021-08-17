#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "Hit2.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Vertex2.h"
#include "Retta.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom.h"
#define RISOLUZIONE 1

double recons(TClonesArray *, TClonesArray* , TH1D*);

void LeggiTree(){
  gRandom ->SetSeed(0);
  // definizione classe vertice
  Vertex2 *ptrvrt = nullptr;
  // Dichiarazione TClonesArray
  // Dichiarazione del TConsArray per gestire gli hit. Per ora ne dichiaro 3.
  // beam pipe
  TClonesArray *Hit_BP = new TClonesArray("Hit2",100);
  TClonesArray& hbp = *Hit_BP;
  // layer 1
  TClonesArray *Hit_L1 = new TClonesArray("Hit2",100);
  TClonesArray& hl1 = *Hit_L1;
  
  // layer 2
  TClonesArray *Hit_L2 = new TClonesArray("Hit2",100);
  TClonesArray& hl2 = *Hit_L2;
  //Apertura file di input
  TFile hfile("htree.root");
  
  //Lettura TTree  e branch
  TTree *tree = (TTree*)hfile.Get("T");
  TBranch *b1=tree->GetBranch("TVert");
  TBranch *b2=tree->GetBranch("BeamPipe");
  TBranch *b3=tree->GetBranch("Layer1");
  TBranch *b4=tree->GetBranch("Layer2");
 




  // Definizione degli indirizzi per la lettura dei dati su ttree
  b1->SetAddress(&ptrvrt);
  b2->SetAddress(&Hit_BP);
  b3->SetAddress(&Hit_L1);
  b4->SetAddress(&Hit_L2);
  

  TH1D* efficienza = new TH1D("histeff", "efficienza", 101, -0.5, 100.5);
  TH1D* residui = new TH1D("histTest", "valori in z", 100, -0.5, 0.5);
  TH1D* tracklet = new TH1D("tracce", "histo", 100, -25., 25.);
  TH1D* molteplicita = new TH1D("histmolt", "molteplicita'", 101, -0.5, 100.5);
  TH1D* evtric = new TH1D("histmolt2", "ev rec'", 101, -0.5, 100.5);
  double RL1 = 4., RL2 = 7.;
  double deltazeta = 0.12;
  double deltarphi = 0.03;
  // loop sugli ingressi nel TTree
  for(int ev = 0; ev< tree -> GetEntries(); ev++){
  tree->GetEvent(ev);
  molteplicita -> Fill(ptrvrt->GetMult());
  }
  new TCanvas("c0","molteplicità",600,600);
  molteplicita -> DrawCopy();

  for(int ev=0;ev<tree->GetEntries();ev++){
    if(10*ev%(tree ->GetEntries()) == 0) cout<<"Sto processando l'evento "<<ev<<endl;
    tree->GetEvent(ev);
    double Zvalue = recons(Hit_L2, Hit_L1, tracklet);
    tracklet -> Reset();
    double residuo = Zvalue - (ptrvrt -> GetZ());
    residui -> Fill(residuo);
    
    if(fabs(residuo) <= 0.75 ) {
    int mult = ptrvrt -> GetMult();
    evtric -> Fill(mult);
    }
  }
  new TCanvas("cuwu","evtric",600,600);
  evtric -> DrawCopy();
  
  for(int i = 0; i < 101; i++){
    int mult = molteplicita -> GetBinContent(i);
    if(mult != 0){
        int eventiric = evtric -> GetBinContent(i);
        efficienza -> Fill(i, double(eventiric)/double(mult));
    }
  }
  new TCanvas("cuwu","evtric",600,600);
  evtric -> DrawCopy();
  
  new TCanvas("c1","residui",600,600);
  residui ->DrawCopy(); 
  new TCanvas("c2","efficienza",600,600);
  efficienza -> DrawCopy();
}


double recons(TClonesArray *Hit_L2, TClonesArray *Hit_L1, TH1D *tracklet){
  double RL1 = 4., RL2 = 7.;
  double deltazeta = 0.12;
  double deltarphi = 0.03;
  TClonesArray& hl1 = *Hit_L1;
  TClonesArray& hl2 = *Hit_L2;

  for(int l =0; l < hl2.GetEntries(); l++){
      Hit2 *tst2 = (Hit2*)Hit_L2->At(l);
      if(tst2 -> GetStatus() != -999){
        double Phi2 = tst2 ->GetPhi();
        double Z2 = tst2 -> GetZP();
        #if RISOLUZIONE
        // SIMULAZIONE DELLA RISOLUZONE DEL RIVELATORE
        Z2 = gRandom->Gaus(Z2,deltazeta);
        Phi2 =  gRandom->Gaus(Phi2,deltarphi/RL2);  
        #endif   
        double counter = tst2 -> GetCounter();
        for(int k = 0; k < hl1.GetEntries(); k++){
          Hit2 *tst1 = (Hit2*)Hit_L1->At(k);
          if((tst1 -> GetStatus() != -999) /*&& (counter == tst1 -> GetCounter())*/){
            double Phi1 = tst1 -> GetPhi();
            double Z1 = tst1 -> GetZP();
            #if RISOLUZIONE
            // SIMULAZIONE DELLA RISOLUZONE DEL RIVELATORE
            Z1 = gRandom->Gaus(Z1,deltazeta);
            Phi1 =  gRandom->Gaus(Phi1,deltarphi/RL1);  
            #endif
            if(fabs(Phi2 - Phi1) < 0.1){
              double m = (RL2 - RL1)/(Z2 - Z1);
              double q = (RL1 - m*Z1);
              double  Zrec = - q / m;
              tracklet -> Fill(Zrec);
            } 
          }   
        }  
      }
    }
    int binmax = tracklet->GetMaximumBin(); 
    double Zvalue = tracklet ->GetXaxis()->GetBinCenter(binmax);
    return Zvalue;
  }