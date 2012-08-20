#define exercise1_cxx
#include "exercise1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#ifdef __CINT__
int exercise1() {
  TFile *f = new TFile("pgs_events_Q1.root");
  TTree *tree = (TTree*)gDirectory->Get("LHCO");
//
.L exercise1.C;
//
//  exercise1 t(tree); t.Loop(); return 0; }
  exercise1 t; t.Loop(); return 0; }
#endif

void exercise1::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L exercise1.C
//      Root > exercise1 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   //nentries = fChain->GetEntriesFast()*.1;		// data for first 10% of entries is to be taken
   //nentries = fChain->GetEntriesFast()*.01;		// data for first 1%  of entries is to be taken

  double minf=-0.5;
  double maxf=4000.5;
  int range=maxf-minf;
  TH1F *lorentzVectorMvalues = new TH1F("lorentzVectorMvalues", "Total mass of 2 electron events",
				range,minf,maxf);

   unsigned int n2eleevents = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      if ( Electron_size != 2 ) continue;

	/*
	Q: Could you modify the code such that only electrons of opposite charge are considered in the reconstruction?
	*/
//      if(Electron_Charge[0]==Electron_Charge[1]) continue;

      n2eleevents++;
      TLorentzVector el1, el2;
      el1.SetPtEtaPhiM( Electron_PT[0], Electron_Eta[0], Electron_Phi[0], 0.0005 ); // in GeV
      el2.SetPtEtaPhiM( Electron_PT[1], Electron_Eta[1], Electron_Phi[1], 0.0005 ); // in GeV

      TLorentzVector res = el1+el2;
      cout << res.M() << "\t";
//      cout << res.M() << "\n";
    
      double tmp=res.M();
      lorentzVectorMvalues->Fill(tmp);
	
   }
   cout << "\nNumber of events with exactly 2 electrons = " << n2eleevents << endl;

////    TF1 *gfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
     TF1 *gfit = new TF1("gfit","gaus",minf,maxf); // Create the fit function
//     lorentzVectorMvalues->Fit("gfit","LL");	//Log-likelihood turu

     /////////////AAAAAAAAAAAAAAAAAAAAAAAAAAAAA/////////////////////////////
   ///// 8.7.2012
   double min1=850;
   double max1=1050;
   TF1 *fit1=new TF1("fit1","gaus",min1,max1); // Create the fit function 
   TF1 *fitpol1=new TF1("fitpol1","pol1",min1,max1); // Create the fit function 
   TF1 *fitpol2=new TF1("fitpol2","pol2",min1,max1); // Create the fit function
//   TF1 *fitgausspol1=new TF1("fitgausspol1","gaus+pol1",min1,max1); // Create the fit function  
//   TF1 *fitgausspol2=new TF1("fitgausspol2","gaus+pol2",min1,max1); // Create the fit function  
//   TF1 *fitgausspol3=new TF1("fitgausspol3","gaus+pol3",min1,max1); // Create the fit function  
//   TF1 *fitgausspol11=new TF1("fitgausspol11","fit1+fitpol1",min1,max1); // Create the fit function  
//   TF1 *fitgausspol22=new TF1("fitgausspol22","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x+[5]*x^2",min1,max1); // Create the fit function  
//   TF1 *fitgausspol111=new TF1("fitgausspol111","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",min1,max1); // Create the fit function

   //////////////////////////////////////////////////////////////////////////////
   
   //////////////////////////////////////////////////////////////////////////////
   // Define the parameter array for the total function
   Double_t par1[5];
   TF1 *fitnoise1gaus=new TF1("fitnoise1gaus","gaus",min1,max1); 
   TF1 *fitnoise1pol=new TF1("fitnoise1pol","pol1",min1,max1);
   TF1 *fitnoise1total=new TF1("fitnoise1total","gaus(0)+pol1(3)",min1,max1);
   
   cout<<endl<<"*****\t"<<fitnoise1total->GetName()<<"\t*****";
   cout << endl << "\t*****\t" << fitnoise1gaus->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   //cout << "-----" << "fitnoise1gaus" << "-----" << endl;
   lorentzVectorMvalues->Fit("fitnoise1gaus","R+");
   cout << endl <<"\t*****\t" << fitnoise1pol->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fitnoise1pol","R+");
   
   fitnoise1gaus->GetParameters(&par1[0]);
   fitnoise1pol->GetParameters(&par1[3]);
   fitnoise1total->SetParameters(&par1[0]);
   
   fitnoise1total->SetLineColor(2);
   cout << endl << "\t*****\t" << fitnoise1total->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fitnoise1total","R+");
   //////////////////////////////////////////////////////////////////////////////
   
   //////////////////////////////////////////////////////////////////////////////
   // Define the parameter array for the total function
      Double_t par2[6];
   TF1 *fitnoise2gaus=new TF1("fitnoise2gaus","gaus",min1,max1);
   TF1 *fitnoise2pol=new TF1("fitnoise2pol","pol2",min1,max1);
   TF1 *fitnoise2total=new TF1("fitnoise2total","gaus(0)+pol2(3)",min1,max1);

   cout<<endl<<"*****\t"<<fitnoise2total->GetName()<<"\t*****";
   cout << endl << "\t*****\t" << fitnoise2gaus->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fitnoise2gaus","R+");
   cout << endl << "\t*****\t" << fitnoise2pol->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fitnoise2pol","R+");
   
   fitnoise2gaus->GetParameters(&par2[0]);
   fitnoise2pol->GetParameters(&par2[3]);
   fitnoise2total->SetParameters(&par2[0]);
   
   fitnoise2total->SetLineColor(3);
   cout << endl << "\t*****\t" << fitnoise2total->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fitnoise2total","R+");
//////////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////
   /*// Define the parameter array for the total function
      Double_t par3[8];
   TF1 *fitnoise3gaus=new TF1("fitnoise3gaus","gaus",min1,max1);
   TF1 *fitnoise3pol=new TF1("fitnoise3pol","pol2",min1,max1);
   TF1 *fitnoise3expo=new TF1("fitnoise3expo","expo",min1,max1);
   TF1 *fitnoise3total=new TF1("fitnoise3total","gaus(0)+pol2(3)+expo(6)",min1,max1);
   
   lorentzVectorMvalues->Fit("fitnoise3gaus","R+");
   lorentzVectorMvalues->Fit("fitnoise3pol","R+");
   lorentzVectorMvalues->Fit("fitnoise3expo","R+");
   
   fitnoise3gaus->GetParameters(&par3[0]);
   fitnoise3pol->GetParameters(&par3[3]);
   fitnoise3expo->GetParameters(&par3[6]);
   fitnoise3total->SetParameters(&par3[0]);
   
   fitnoise3total->SetLineColor(4);
   lorentzVectorMvalues->Fit("fitnoise3total","R+");*/
//////////////////////////////////////////////////////////////////////////////   
   
   
   fit1->SetLineColor(1);
   fitpol1->SetLineColor(6);
   fitpol2->SetLineColor(7);
   //   fitgausspol1->SetLineColor(6);
   //   fitgausspol2->SetLineColor(7);
   //   fitgausspol3->SetLineColor(8);
   //fitgausspol11->SetLineColor(9);
//   fitgausspol111->SetLineColor(5);
//   fitgausspol22->SetLineColor(9);
   fitnoise1total->SetLineColor(2);
   fitnoise2total->SetLineColor(3);
   /*
   1	black
   2	red
   3	light green
   4	blue
   5	yellow
   6	magenta
   7	cyan
   8	green*/
   
   cout << endl << "\t*****\t" << fit1->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fit1","R+");	
   
   cout << endl << "\t*****\t" << fitpol1->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fitpol1","R+");
   
   cout << endl << "\t*****\t" << fitpol2->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
   lorentzVectorMvalues->Fit("fitpol2","R+");	
//   lorentzVectorMvalues->Fit("fitgausspol1","R+");
//   lorentzVectorMvalues->Fit("fitgausspol2","R+");
//   lorentzVectorMvalues->Fit("fitgausspol3","R+");
//   lorentzVectorMvalues->Fit("fitgausspol11","R+");
//   lorentzVectorMvalues->Fit("fitgausspol111","R+");
//   lorentzVectorMvalues->Fit("fitgausspol22","R+");
   
   // legend ekle.
   leg_hist = new TLegend(0.7,0.4,0.95,0.6);
   
   // add the values of parameters to the legend
   
   TString strGetParams="with ";		// equivalent : TString strGetParams("with");
   for (int i=0;i<3;i++)
   {
	   strGetParams+="p";
	   strGetParams+=i;
	   strGetParams+=" = ";
	   //strGetParams+=par1[i];
	   strGetParams+=+Form("%.3f",par1[i]);		// 3 numbers after the decimal point (virgülden sonra 3 basamak)
	   strGetParams+=" , ";
	   //strGetParams+="p"+i+" = ";		// must concatenate one by one ( böyle bir kod çalýþmýyor, teker teker eklemek gerekiyor.) 
	   //strGetParams<<"p"<<i+" = ";//"+par1[i]+" , ";  // must concatenate one by one ( böyle bir kod çalýþmýyor, teker teker eklemek gerekiyor.)
   }
   
   leg_hist->AddEntry(fit1,"gaus = p0*exp(-0.5*((x-p1)/p2)^2)","l");
   leg_hist->AddEntry(fitpol1,"pol1 = p0 + p1*x","l");
   leg_hist->AddEntry(fitpol2,"pol2 = p0 + p1*x + p2*x^2","l");
   //leg_hist->AddEntry(fitnoise1total,"gaus+pol1 "+strGetParams,"l");	// given up, makes the legend unnecessarily wide
   leg_hist->AddEntry(fitnoise1total,"gaus+pol1 ","l");
   leg_hist->AddEntry(fitnoise2total,"gaus+pol2","l");
   leg_hist->Draw();
   
   
   /////////////AAAAAAAAAAAAAAAAAAAAAAAAAAAAA/////////////////////////////
   
    TCanvas *c1 = new TCanvas();
    c1->Divide(2);

   c1->cd(1);
   lorentzVectorMvalues->Draw();

    c1->cd(2);    
    LHCO->Draw("Electron.PT");
    
    // ZOOM IN AREA OF INTEREST
    TH1F *lorentzVectorMvaluesSubRange=lorentzVectorMvalues->Clone();
    lorentzVectorMvaluesSubRange->SetName("lorentzVectorMvaluesSubRange");
    lorentzVectorMvaluesSubRange->SetTitle(" with range of interest");
    lorentzVectorMvaluesSubRange->GetXaxis()->SetRange(min1,max1);
     
    TCanvas *c2 = new TCanvas();
    lorentzVectorMvaluesSubRange->Draw();
    leg_hist->Draw();		// "legend" etiketlerini çiz.
}

/*
   unsigned int n2eleevents = 0;
...
      if ( Electron_size != 2 ) continue;
      n2eleevents++;
...
   cout << "\nNumber of events with exactly 2 electrons = " << n2eleevents << endl;
*/

/*
      TLorentzVector el1, el2;
      el1.SetPtEtaPhiM( Electron_PT[0], Electron_Eta[0], Electron_Phi[0], 0.0005 ); // in GeV
      el2.SetPtEtaPhiM( Electron_PT[1], Electron_Eta[1], Electron_Phi[1], 0.0005 ); // in GeV

      TLorentzVector res = el1+el2;
      cout << res.M() << "\t";
*/
