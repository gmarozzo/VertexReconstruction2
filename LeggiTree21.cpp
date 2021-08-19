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
#include "TVector.h"


#define RISOLUZIONE 1

double recons(TClonesArray *, TClonesArray* , TH1D*);
double recons2(TClonesArray *, TClonesArray*, double &varianza);

void LeggiTree(){
  gErrorIgnoreLevel = kError;
  vector <double> Zricostruiti;
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

  TH1D* tracklet = new TH1D("tracce", "histo", 100, -25., 25.);
  TH1D* residui = new TH1D("histTest", "valori in z", 100, -0.15, 0.15); 

  TH1D* efficienza = new TH1D("histeff", "efficienza vs mult", 100, 0.5, 100.5); 
  TH1D* risoluzione = new TH1D("hrisol", "risoluzione vs mult", 100, 0.5, 100.5);
  TH1D* molteplicita = new TH1D("histmolt", "molteplicita'", 100, 0.5, 100.5);
  TH1D* evtric = new TH1D("histmolt2", "ev rec'", 100, 0.5, 100.5);

  TH1D* efficienzaz = new TH1D("histeff2", "efficienza vs z vertice", 10, -20, 20);
  TH1D* risoluzionez = new TH1D("hrisol2", "risoluzione vs z vertice", 10, -20, 20);
  TH1D* zvertice = new TH1D("hz", "z vertice generato", 10, -20, 20);
  TH1D* zvertric = new TH1D("hzric", "z vertice ricostruito", 10, -20, 20);
  
  double RL1 = 4., RL2 = 7.;
  double deltazeta = 0.012;
  double deltarphi = 0.003;
  float sommaresidui[100];
  for(int i = 0; i < 100; i++){
      sommaresidui[i] = 0.;
  }
  /*// loop sugli ingressi nel TTree
  for(int ev = 0; ev< tree -> GetEntries(); ev++){
  tree->GetEvent(ev);
  molteplicita -> Fill(ptrvrt->GetMult());
  }
  */

  for(int ev=0;ev<tree->GetEntries();ev++){
    double varianza = 0.;
    if(10*ev%(tree ->GetEntries()) == 0) cout<<"Sto processando l'evento "<<ev<<endl;
    tree->GetEvent(ev);
    double Zvalue = recons2(Hit_L2, Hit_L1, varianza);
    //tracklet -> Reset("M");
    float residuo = Zvalue - (ptrvrt -> GetZ());
    residui -> Fill(residuo);
    int mult = ptrvrt -> GetMult();
    double zvert = ptrvrt -> GetZ();
    molteplicita -> Fill(mult);
    zvertice-> Fill(zvert);
    if(fabs(residuo) <= 3*sqrt(varianza)){
    evtric -> Fill(mult);
    zvertric -> Fill(zvert);
    if(mult >= 1) sommaresidui[mult - 1] += (residuo*residuo);
    if(fabs(sommaresidui[mult - 1]) > 1e5) cout<<"!!WARNING!!"<<residuo<<" , "<<ptrvrt -> GetMult()<<endl;
    risoluzionez->Fill(zvert,residuo*residuo);
    }
  }
  new TCanvas("c0","evtric",600,600);
  //evtric -> DrawCopy();
  zvertric -> DrawCopy();
  
  new TCanvas("c1","molteplicità",600,600);
  //molteplicita -> DrawCopy();
  zvertice -> DrawCopy();

  for(int i = 0; i < 100; i++){
    int mult = molteplicita -> GetBinContent(i + 1);
    if(mult != 0){
    int eventiric = evtric -> GetBinContent(i + 1);
    efficienza -> SetBinContent(i + 1, double(eventiric)/double(mult));}
  }

  for(int i = 0; i < 10; i++){
    double zvert = zvertice -> GetBinContent(i+1);
    if(zvert != 0){
    double zverticeric = zvertric -> GetBinContent(i+1);
    efficienzaz -> SetBinContent(i+1,(double)zverticeric/zvert);
    }
   }
    
  
  
  for(int i = 0; i < 100; i++){
    int eventi = evtric -> GetBinContent(i + 1);
    if(eventi > 1e-12){
      risoluzione -> Fill(i + 1, sqrt(sommaresidui[i]/(eventi-1)));
    //cout<<"risoluzione "<< sqrt(sommaresidui[i]/eventi)<<endl;
    //cout<<"sommaresidui"<< sommaresidui[i]<<endl;
    }
  }

  for(int i = 0; i < 10; i++){
    int eventi = zvertric -> GetBinContent(i + 1);
    if(eventi > 1){
      risoluzionez -> SetBinContent(i + 1, sqrt(risoluzionez->GetBinContent(i+1)/(eventi-1)));
      cout<<risoluzionez->GetBinContent(i+1)<<endl;
    }
  }
  


  new TCanvas("c2","residui",600,600);
  residui -> DrawCopy(); 
  
  new TCanvas("c3","efficienza",600,600);
  efficienza -> SetMarkerStyle(20);
  efficienza -> SetMarkerSize(1);
  efficienza -> DrawCopy("histp");
  
  new TCanvas("c5","efficienzaz",600,600);
  efficienzaz -> SetMarkerStyle(20);
  efficienzaz -> SetMarkerSize(1);
  efficienzaz -> DrawCopy("histp");
  new TCanvas("c4","risoluzione",600,600);
  risoluzione -> SetMarkerStyle(20);
  risoluzione -> SetMarkerSize(1);
  risoluzione -> DrawCopy("histp");
  new TCanvas("c6","risoluzione",600,600);
  risoluzionez -> SetMarkerStyle(20);
  risoluzionez -> SetMarkerSize(1);
  risoluzionez -> DrawCopy("histp");
}


double recons(TClonesArray *Hit_L2, TClonesArray *Hit_L1, TH1D *tracklet){
  double RL1 = 4., RL2 = 7.;
  double deltazeta = 0.012;
  double deltarphi = 0.003;
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
            if(fabs(Phi2 - Phi1) < 1.){
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

  double recons2(TClonesArray *Hit_L2, TClonesArray *Hit_L1, double &varianza){
  double RL1 = 4., RL2 = 7.;
  double deltazeta = 0.012;
  double deltarphi = 0.003;
  TClonesArray& hl1 = *Hit_L1;
  TClonesArray& hl2 = *Hit_L2;
  TH1D* tracklet = new TH1D("tracce", "histo", 100, -20., 20.);
  vector <double> Zricostruiti;
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
              Zricostruiti.push_back(Zrec);
            } 
          }   
        }  
      }
    }
    int binmax = tracklet->GetMaximumBin(); 
    double Zvalue = tracklet ->GetXaxis()->GetBinCenter(binmax);
    double ZvalueP = tracklet -> GetXaxis()->GetBinCenter(binmax + 3);
    double ZvalueM = tracklet -> GetXaxis()->GetBinCenter(binmax - 3);
    tracklet -> GetXaxis() -> SetRange(binmax-3, binmax+3);
    double Media = tracklet -> GetMean();
    double MediaVector = 0.;
    int contatore = 0;
    //BRUTTO DA RIVEDERE
    for(int i = 0; i < Zricostruiti.size(); i++){
      if (fabs(Zricostruiti[i] - Media) <= 0.25){
      MediaVector += Zricostruiti[i];
      contatore ++;
      }
    }
    MediaVector = MediaVector/contatore;
    
    for(int i = 0; i < Zricostruiti.size(); i++){
      if (fabs(Zricostruiti[i] - Media) <= 0.25){
      double x = (MediaVector - Zricostruiti[i]);
      varianza += x*x/(contatore - 1);
      }
    }
    
    //tracklet -> GetXaxis() -> SetRange(0, 100);
    return MediaVector;
  }
