#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "Vertex2.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "MyRandom3.h"
#include "Retta.h"
#include "Hit2.h"
#include "TH1.h"
#include "TCanvas.h"
#include <TH3D.h>


#define DEBUG 0 // stampa a schermo degli eventi (da implementare)
#define DEBUG2 0 // controllo degli eventi (istogramma 3-dim)
#define DEBUG3 0 // controllo degli eventi fuori dal range
#define DEBUG4  0 // controllo degli eventi di fondo
#define DEBUG5 0 // controllo delle tracce (bug nel debug >-<)
#define DEBUG6 0 // prova per ricostruzione tracce
#define SCTRM 1 // scattering multiplo
#define RISOLUZIONE 1 // applica gli effetti di risoluzione del rivelatore


TH1F* maniphist();
void ROTATE(Double_t , Double_t , Double_t , Double_t , Double_t *cd);
void multiplescattering(double *angs, double *ang, double *cd);

void GeneraTree(){
  
  TH1D* residui = new TH1D("histTest", "valori in z", 50, -0.1, 0.1);
  TH1D* zetaveri = new TH1D("hhh","zeta generati",100,-5,5);
  TH3D* palla = new TH3D("palla","sfera",100,-137.1,137.1,100,-137.1,137.1,100,-137.1,137.1);
  TH1D* tracket = new TH1D("tracce", "histo", 100, -5., 5.);
  double X, Y, Z;
  double RBP = 3.;
  double RL1 = 4.;
  double RL2 = 7.;
  double Zmin1 = -13.5, Zmax1 = 13.5;
  double Zmin2 = -13.5, Zmax2 = 13.5;
  double deltazeta = 0.12;
  double deltarphi = 0.03;
  int nfondo = 10;
  int nev = 1e4;
  
  TFile *f1 = new TFile("kinem.root");
  TH1F *pseta = (TH1F*) f1 -> Get("heta");
  pseta -> SetDirectory(0);
  f1 -> Close();
  
  TH1F *heta = maniphist();
  //heta -> Draw();
  
  // Apertura del file di output
  TFile hfile("htree.root","RECREATE");
  // dichiarazione del TTree
  TTree *tree = new TTree("T","TTree con 4 branches");
  // se invertissi l'ordine dovrei scrivere
  // tree->SetDirectory(&hfile);

  // Dichiarazione del vertice. Puntatore nullo ad un oggetto vertice
  Vertex2* ptrvrt = nullptr;
  
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
  


  // Dichiarazione del branch del TTree
  tree -> Branch("TVert", &ptrvrt);
  tree -> Branch("BeamPipe", &Hit_BP);
  tree -> Branch("Layer1", &Hit_L1);
  tree -> Branch("Layer2", &Hit_L2);
   
  
  
  
  
  //INIZIA IL LOOP SUGLI EVENTI
  for(int i=0; i<nev;i++){ // loop sugli eventi
    if (i%(nev/10) == 0)cout<<"Sto processando l' "<< i<<"  esimo evento"<<endl;
    // Genero una molteplicità e un vertice
    int numpart=0;
    //while(numpart<=0){
      numpart=(int)(0.5+gRandom->Gaus(50.,20.));
    //}
    
    numpart = int((gRandom -> Rndm())*100);
    X=gRandom->Gaus(0,0.01);
    Y=gRandom->Gaus(0,0.01);
    Z=gRandom->Gaus(0.,5.3);
    zetaveri->Fill(Z);
    
    Vertex2 *vrt = new  Vertex2(X, Y, Z, numpart);
    ptrvrt = vrt;  
    
    // Generato il vertice passo alla macro che gestisce gli hit.
    // Un hit viene creato a partire da un angolo theta, da un angolo phi e da un counter + posizione x, y, z
    for(int j = 0; j < numpart; j++){
      double X0_L1[3];
      double X0_L2[3];
      double ang[2];
      double position[3];
      //ANGOLI DI PARTENZA
      double theta = pseta -> GetRandom();
      double phi = (gRandom -> Rndm())*2.*M_PI;
      ang[0]=theta;
      ang[1]=phi;
      position[0]=X;
      position[1]=Y;
      position[2]=Z;
      Retta *retta1 = new Retta(position,ang);
      retta1->intpoint(position, RBP);   
      new(hbp[j])Hit2(ang[0], ang[1], position, j);
      //Hit2 *tst0=(Hit2*)Hit_BP->At(j);
      //if (fabs(position[2]) >= 27) tst0 ->SetStatus(-999);
      
      // LAYER 1
      X0_L1[0] = position[0];
      X0_L1[1] = position[1];
      X0_L1[2] = position[2];
      if(SCTRM){
      retta1->multiplescattering();
      ang[0]=retta1->GetTheta();
      ang[1]=retta1->GetPhi();
      }
      Retta *retta2 = new Retta(position,ang);
      retta2->intpoint(position, RL1); 
      new(hl1[j])Hit2(ang[0], ang[1], position, j );
      Hit2 *tst1=(Hit2*)Hit_L1->At(j);
      if (fabs(position[2]) >= 13.5) tst1 ->SetStatus(-999);
      #if RISOLUZIONE
      // SIMULAZIONE DELLA RISOLUZONE DEL RIVELATORE
      double zetars = gRandom->Gaus(tst1 -> GetZP(),deltazeta);
      double phirs =  gRandom->Gaus(tst1 -> GetPhi(),deltarphi/RL1);
      tst1 -> Rec_Z = zetars;
      tst1 -> Rec_Phi = phirs; 
      #endif
      // LAYER 2
      X0_L2[0] = position[0];
      X0_L2[1] = position[1];
      X0_L2[2] = position[2];
      if(SCTRM){
      retta2->multiplescattering();
      ang[0]=retta2->GetTheta();
      ang[1]=retta2->GetPhi();
      }
      Retta *retta3 = new Retta(position,ang);
      retta3->intpoint(position, RL2); 
      new(hl2[j])Hit2(ang[0], ang[1], position, j);
      Hit2 *tst2=(Hit2*)Hit_L2->At(j);
      if (fabs(position[2]) >= 13.5) tst2 ->SetStatus(-999);
      #if RISOLUZIONE
      // SIMULAZIONE DELLA RISOLUZIONE DEL RIVELATORE
      zetars = gRandom->Gaus(tst2 -> GetZP(),deltazeta);
      phirs =  gRandom->Gaus(tst2 -> GetPhi(),deltarphi/RL2);
      tst2 -> Rec_Z = zetars;
      tst2 -> Rec_Phi = phirs;  
      #endif   
    }


    // SIMULAZIONE DI HIT CASUALI
    for(int i = numpart; i < numpart + nfondo; i++){
    double ang[2]; // angolo dell'hit rispetto al lab
    double position[3]; // Posizione degli hit

    // simulazione del fondo sul layer 1
    ang[1] = (gRandom -> Rndm())*2.*M_PI;
    position[2] = Zmin1 + (gRandom -> Rndm())*(Zmax1 - Zmin1);
    
    ang[0] = TMath::ATan(RL1/position[2]);
    position[0] = RL1 * TMath::Sin(ang[0])*TMath::Cos(ang[1]);
    position[1] = RL1 * TMath::Sin(ang[0])*TMath::Sin(ang[1]);
    
    new(hl1[i])Hit2(ang[0], ang[1], position, -999);
    
    /*
    // SIMULAZIONE DELLA RISOLUZIONE DEL RIVELATORE
    Hit2 *tst1=(Hit2*)Hit_L1->At(i);
    tst1 -> Rec_Z = position[2];
    tst1 -> Rec_Phi = ang[1];
    */
    // simulazione del fondo sul layer 2
    ang[1] = (gRandom -> Rndm())*2.*M_PI;
    position[2] = Zmin2 + (gRandom -> Rndm())*(Zmax2 - Zmin2);
    
    ang[0] = TMath::ATan(RL2/position[2]);
    position[0] = RL2 * TMath::Sin(ang[0])*TMath::Cos(ang[1]);
    position[1] = RL2 * TMath::Sin(ang[0])*TMath::Sin(ang[1]);
    
    new(hl2[i])Hit2(ang[0], ang[1], position, -999);
    /*
    // SIMULAZIONE DELLA RISOLUZONE DEL RIVELATORE
    Hit2 *tst2=(Hit2*)Hit_L2->At(i);
    tst2 -> Rec_Z = position[2];
    tst2 -> Rec_Phi = ang[1];
    */
    }
    
    
    #if DEBUG
    // Debug
    
    if(i%10 == 0){
    char titolo[50];
    char canvas[50];
    sprintf(canvas, "cv %d", i/10 + 1);
    sprintf(titolo, "evento %d", i);
    double x0 = ptrvrt -> GetX();
    double y0 = ptrvrt -> GetY();
    double z0 = ptrvrt -> GetZ();
    for (int j=0; j<hbp.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_BP->At(j);
      
      double th = tst -> GetTheta();
      double pphi = tst -> GetPhi();
      double c[3];
      c[0]=TMath::Sin(th)*TMath::Cos(pphi);
      c[1]=TMath::Sin(th)*TMath::Sin(pphi);
      c[2] = TMath::Cos(th);
      double t = intersection(x0, y0, c, RBP);
      palla->Fill(x0 + t*c[0],y0 + t*c[1], z0 + t*c[2]);
    }
    new TCanvas(canvas,titolo,600,600);
    palla -> DrawCopy();
    palla -> Reset();
    }
    // fine del debug
    #endif
    
    #if DEBUG2
    // Debug
    
    if(i%10 == 0){
    char titolo[50];
    char canvas[50];
    sprintf(canvas, "cv %d", i/10 + 1);
    sprintf(titolo, "evento %d", i);
    for (int j=0; j<hbp.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_BP->At(j);
      double xx = tst -> GetXP();
      double yy = tst -> GetYP();
      double zz = tst -> GetZP();
      palla->Fill(xx, yy, zz);
    }
    for (int j=0; j<hl1.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_L1->At(j);
      double xx = tst -> GetXP();
      double yy = tst -> GetYP();
      double zz = tst -> GetZP();
      palla->Fill(xx, yy, zz);
    }
    for (int j=0; j<hl2.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_L2->At(j);
      double xx = tst -> GetXP();
      double yy = tst -> GetYP();
      double zz = tst -> GetZP();
      palla->Fill(xx, yy, zz);
      }
    
    new TCanvas(canvas,titolo,600,600);
    palla -> SetMarkerSize(0.5);
    palla -> SetMarkerStyle(20);
    palla -> DrawCopy();
    palla -> Reset();
    }
    // fine del debug
    #endif
    
    #if DEBUG3
    //debug
    for (int j=0; j<hbp.GetEntries(); j++){
      Hit2 *tst001=(Hit2*)Hit_BP->At(j);
      if(tst001 -> GetStatus() == -999) residui -> Fill(tst001 -> GetZP());
    }
    //fine del debug
    #endif

    #if DEBUG4
    // Debug
    
    if(i%10 == 0){
    char titolo[50];
    char canvas[50];
    sprintf(canvas, "cv %d", i/10 + 1);
    sprintf(titolo, "evento %d", i);
    
    /*
    for (int j=0; j<hbp.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_BP->At(j);
      if (tst -> GetCounter() == -999){
        double xx = tst -> GetXP();
        double yy = tst -> GetYP();
        double zz = tst -> GetZP();
        palla->Fill(xx, yy, zz);
      }
      
    }
    */
    for (int j=0; j<hl1.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_L1->At(j);
      if (tst -> GetCounter() == -999){
        double xx = tst -> GetXP();
        double yy = tst -> GetYP();
        double zz = tst -> GetZP();
        palla->Fill(xx, yy, zz);
      }
    }
    for (int j=0; j<hl2.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_L2->At(j);
      if (tst -> GetCounter() == -999){
        double xx = tst -> GetXP();
        double yy = tst -> GetYP();
        double zz = tst -> GetZP();
        palla->Fill(xx, yy, zz);
      }
    }
    /*
    new TCanvas(canvas,titolo,600,600);
    palla -> SetMarkerSize(0.5);
    palla -> SetMarkerStyle(20);
    palla -> DrawCopy();
    */
    palla -> Reset();
    
    }
    // fine del debug
    #endif
    
     #if DEBUG5
    // Debug
    
    if(i%10 == 0){
    char titolo[50];
    char canvas[50];
    sprintf(canvas, "cv %d", i/10 + 1);
    sprintf(titolo, "evento %d", i);
    new TCanvas(canvas,titolo,600,600);
    palla -> SetMarkerSize(0.5);
    palla -> SetMarkerStyle(20);
    for (int j=0; j<hbp.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_BP->At(j);
      double xx = tst -> GetXP();
      double yy = tst -> GetYP();
      double zz = tst -> GetZP();
      palla->Fill(xx, yy, zz);
      if(tst -> GetCounter() == 0) palla -> SetMarkerColor(1);
      if(tst -> GetCounter() == 1) palla -> SetMarkerColor(2);
      if(tst -> GetCounter() == 2) palla -> SetMarkerColor(3);
      palla -> SetMarkerSize(0.5);
      palla -> SetMarkerStyle(20);
      palla -> DrawCopy("Same");
      palla -> Reset(); 
    }
    for (int j=0; j<hl1.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_L1->At(j);
      double xx = tst -> GetXP();
      double yy = tst -> GetYP();
      double zz = tst -> GetZP();
      palla->Fill(xx, yy, zz);
      if(tst -> GetCounter() == 0) palla -> SetMarkerColor(1);
      if(tst -> GetCounter() == 1) palla -> SetMarkerColor(2);
      if(tst -> GetCounter() == 2) palla -> SetMarkerColor(3);
      palla -> SetMarkerSize(0.5);
      palla -> SetMarkerStyle(20);
      palla -> DrawCopy("Same");
      palla -> Reset();
    }
    for (int j=0; j<hl2.GetEntries(); j++){
      Hit2 *tst=(Hit2*)Hit_L2->At(j);
      double xx = tst -> GetXP();
      double yy = tst -> GetYP();
      double zz = tst -> GetZP();
      palla->Fill(xx, yy, zz);
      if(tst -> GetCounter() == 0) palla -> SetMarkerColor(1);
      if(tst -> GetCounter() == 1) palla -> SetMarkerColor(2);
      if(tst -> GetCounter() == 2) palla -> SetMarkerColor(3);
      palla -> SetMarkerSize(0.5);
      palla -> SetMarkerStyle(20);
      palla -> DrawCopy("Same");
      palla -> Reset();
      }
    /*new TCanvas(canvas,titolo,600,600);
    palla -> SetMarkerSize(0.5);
    palla -> SetMarkerStyle(20);
    palla -> DrawCopy("Same");
    palla -> Reset();
    */
    }
    // fine del debug
    #endif

    #if DEBUG6
    
    for(int l =0; l < hl2.GetEntries(); l++){
      Hit2 *tst2 = (Hit2*)Hit_L2->At(l);
      if(tst2 -> GetStatus() != -999){
        double Phi2 = tst2 ->GetPhi();
        double Z2 = tst2 -> GetZP();
        double counter = tst2 -> GetCounter();
        for(int k = 0; k < hl1.GetEntries(); k++){
          Hit2 *tst1 = (Hit2*)Hit_L1->At(k);
          if((tst1 -> GetStatus() != -999) /*&& (counter = tst1 -> GetCounter())*/){
            double Phi1 = tst1 -> GetPhi();
            double Z1 = tst1 -> GetZP();
            if(fabs(Phi1-Phi2)<1.){
              double m = (RL2 - RL1)/(Z2 - Z1);
              double q = RL1-m*Z1;
              double  Zrec = - q / m;
              tracket -> Fill(Zrec);
            } 
          }   
        }  
      }
    }
    int binmax = tracket->GetMaximumBin(); 
    double Zvalue = tracket ->GetXaxis()->GetBinCenter(binmax);
    /*if(i%(nev/10) == 0) {
      char canvas[50];
      char titolo[50];
      sprintf(canvas, "cv %d", i/10 + 1);
      sprintf(titolo, "evento %d", i);
      new TCanvas(canvas,titolo,600,600);
      tracket-> DrawCopy();
      cout<<Zvalue<<"  "<<i<<endl;
      cout<<ptrvrt -> GetZ()<<endl;
    }
    */

    
    tracket -> Reset();
    residui -> Fill(Zvalue - (ptrvrt -> GetZ()));
    
    #endif
    tree->Fill();
    Hit_BP->Clear();
    Hit_L1 ->Clear();
    Hit_L2 ->Clear();
  }
  
  new TCanvas("c1","zeta veri",600,600);
  zetaveri ->DrawCopy();
  
  new TCanvas("c2","residui",600,600);
  residui ->DrawCopy();
  // Save all objects in this file
  hfile.Write();

  // Close the file. 
  hfile.Close();
}

TH1F* maniphist(){
  TFile F("kinem.root");
  TH1F *disteta = (TH1F*)F.Get("heta");
  disteta->SetDirectory(0);
  disteta->SetMinimum(0);
  F.Close();
  TAxis *xa=disteta->GetXaxis();
  Double_t step = xa->GetBinWidth(1);
  Int_t b1=xa->FindBin(-2.);
  Int_t b2=xa->FindBin(2.);
  Double_t xlow=xa->GetBinLowEdge(b1);
  Double_t xhig=xa->GetBinUpEdge(b2);
  Int_t nobins=b2-b1+1;
  Double_t step2 = (xhig-xlow)/nobins;
  //cout << "Check: "<<step<<"; "<<step2<<endl;
  TH1F* heta2 = new TH1F("heta2","#eta distribution 2",nobins,xlow,xhig);
  Int_t j=1;
  for(Int_t i=b1;i<=b2;i++)heta2->SetBinContent(j++,disteta->GetBinContent(i));
  //  heta2->Draw();
  //new TCanvas();
  //disteta->Draw();
  //heta2->SetLineColor(2);
  //heta2->Draw("same");
  return heta2;
}