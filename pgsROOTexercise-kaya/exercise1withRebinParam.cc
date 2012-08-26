#define exercise1_cxx
#include "exercise1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#ifdef __CINT__
int exercise1withRebinParam(int rebinParamforRebinFnc) {	//
  TFile *f = new TFile("pgs_events_Q1.root");
  TTree *tree = (TTree*)gDirectory->Get("LHCO");
//
.L exercise1.C;
//
//  exercise1 t(tree); t.Loop(); return 0; }
exercise1 t; t.Loop(rebinParamforRebinFnc); return 0; }
#endif

void exercise1::Loop(int rebinParam)
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

   Long64_t nentriesAll = fChain->GetEntriesFast();
   double percentageOfEvents=.01;
   nentries = fChain->GetEntriesFast()*percentageOfEvents;		// data for first 'percentageOfEvents'*100 of entries is to be taken

  double minf=-0.5;
  double maxf=4000.5;
  int range=maxf-minf;
  TString titleOfTH1F="Total mass of 2 electron events (";
  titleOfTH1F+="first ";
  titleOfTH1F+=Form("%.1f",percentageOfEvents*100);
  titleOfTH1F+="% of ";
  titleOfTH1F+=nentriesAll;
  titleOfTH1F+="entries is studied)";
//  TH1F *lorentzVectorMvalues = new TH1F("lorentzVectorMvalues", "Total mass of 2 electron events",
//				range,minf,maxf);
    TH1F *lorentzVectorMvalues = new TH1F("lorentzVectorMvalues", titleOfTH1F,
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
//      cout << res.M() << "\t";	ekrana yazd�rmak "ASUS netbook" ta �ok zaman al�yor. �imdilik ekrana yazd�rma.
//      cout << res.M() << "\n";
    
      double tmp=res.M();
      lorentzVectorMvalues->Fill(tmp);
   }

   if(rebinParam !=NULL && rebinParam>1)
   {
	   lorentzVectorMvalues->Rebin(rebinParam);		// rebinParam = number of bins to be merged when "Rebin()" is called.
   }
   int numOfBins=lorentzVectorMvalues->GetNbinsX();
   
   
   cout << "\nNumber of events with exactly 2 electrons = " << n2eleevents << endl;
   cout << "\nNumber of bins in 'lorentzVectorMvalues' = " << numOfBins << endl;
   
   TCanvas *c11=new TCanvas();
   
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
   
   //////////////////////////////////////////////////////////////////////////////
   // Define the parameter array for the total function
   Double_t par1[5];
//   TF1 *fitnoise1gaus=new TF1("fitnoise1gaus","gaus",min1,max1); 
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
   
   
   fit1->SetLineColor(1);
   fitpol1->SetLineColor(6);
   fitpol2->SetLineColor(7);
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
	   strGetParams+=+Form("%.3f",par1[i]);		// 3 numbers after the decimal point (virg�lden sonra 3 basamak)
	   strGetParams+=" , ";
	   //strGetParams+="p"+i+" = ";		// must concatenate one by one ( b�yle bir kod �al��m�yor, teker teker eklemek gerekiyor.) 
	   //strGetParams<<"p"<<i+" = ";//"+par1[i]+" , ";  // must concatenate one by one ( b�yle bir kod �al��m�yor, teker teker eklemek gerekiyor.)
   }
   
   leg_hist->AddEntry(fit1,"gaus = p0*exp(-0.5*((x-p1)/p2)^2)","l");
   leg_hist->AddEntry(fitpol1,"pol1 = p0 + p1*x","l");
   leg_hist->AddEntry(fitpol2,"pol2 = p0 + p1*x + p2*x^2","l");
   //leg_hist->AddEntry(fitnoise1total,"gaus+pol1 "+strGetParams,"l");	// given up, makes the legend unnecessarily wide
   leg_hist->AddEntry(fitnoise1total,"gaus+pol1 ","l");
   leg_hist->AddEntry(fitnoise2total,"gaus+pol2","l");
   leg_hist->Draw();
      
   lorentzVectorMvalues->Draw();
   //lorentzVectorMvalues->DrawClone();
   
   /////////////AAAAAAAAAAAAAAAAAAAAAAAAAAAAA/////////////////////////////
   
    TCanvas *c1 = new TCanvas();
    c1->Divide(2);

   c1->cd(1);
   lorentzVectorMvalues->Draw();
/*   When a displayed
   histogram is deleted, its image is automatically removed from the pad. To create a copy of the
   histogram when drawing it, you can use TH1::DrawClone(). This will clone the histogram
   and allow you to change and delete the original one without affecting the clone. */
//   lorentzVectorMvalues->DrawClone();

    c1->cd(2);    
    LHCO->Draw("Electron.PT");
    
    // ZOOM IN AREA OF INTEREST
    TH1F *lorentzVectorMvaluesSubRange=lorentzVectorMvalues->Clone();
        
    lorentzVectorMvaluesSubRange->SetName("lorentzVectorMvaluesSubRange");
    lorentzVectorMvaluesSubRange->SetTitle(" with range of interest");
    //lorentzVectorMvaluesSubRange->GetXaxis()->SetRange(min1,max1);
    lorentzVectorMvaluesSubRange->GetXaxis()->SetRange(min1/rebinParam,max1/rebinParam);	// if "Rebin()" is called, then new range must be specified accordingly.
    
    TCanvas *c2 = new TCanvas();
    
    lorentzVectorMvaluesSubRange->Draw();
    //lorentzVectorMvaluesSubRange->DrawClone();
    leg_hist->Draw();		// "legend" etiketlerini �iz.
    
    //	A PaveStats is a PaveText to draw histogram statistics and fit parameters.
    //	http://root.cern.ch/root/html/TPaveStats.html#TopOfPage
    gPad->Update();
    TPaveStats *st = (TPaveStats*)lorentzVectorMvaluesSubRange->FindObject("stats");
    st->SetOptStat(2211);
    /*
     But in a script file the painting should be forced using gPad->Update() in order to make sure the statistics box is created:

      h->Draw();
      gPad->Update();
      TPaveStats *st = (TPaveStats*)h->FindObject("stats");

Without gPad->Update() the line h->FindObject("stats") returns a null pointer.
To change the type of information for an histogram with an existing TPaveStats one should do:

      st->SetOptStat(mode);
    * */
            
    //Float_t x1a=0.0, y1a=0.9, x2a=0.4, y2a=1.0;
    Double_t x1a=st->GetX1NDC(), y1a=0.0, x2a=st->GetX2NDC(), y2a=st->GetY1NDC();
    y1a=y2a-.1;
    TPaveText *pt2 = new TPaveText(x1a,y1a,x2a,y2a,"brNDC");
    TString ptStr="# bins = \t";    ptStr+=numOfBins;    pt2->AddText(ptStr);
    ptStr="rebin Parameter = \t";	ptStr+=rebinParam;    pt2->AddText(ptStr);
    ptStr="% of events = \t";	ptStr+=Form("%.1f",percentageOfEvents*100);    pt2->AddText(ptStr);
    pt2->Draw();
    

    TFile dosya("different_lorentzVectorMvalues_histos.root","UPDATE");
    TString asd="lvm";
    asd+=rebinParam;
    lorentzVectorMvalues->SetName(asd);
//    lorentzVectorMvalues->SetName("lvm"+rebinParam);
//    lorentzVectorMvalues->SetNameTitle("lvm"+rebinParam);
    lorentzVectorMvalues->Write();
    lorentzVectorMvaluesSubRange->Write();
    dosya.Close();
}


/*
 * What to put inside pavetext
 * 
 *  percentage of events
 *  current	number 0f bins
 *  rebin parameter
 * 
 * */