#define exercise1_cxx
#include "exercise1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#ifdef __CINT__
int exercise1withLL(int rebinParamforRebinFnc) {	//
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

	// VARIABLES TO BE USED IN ADVANCE
	TPaveStats *st;			//	http://root.cern.ch/root/html/TPaveStats.html
	TPaveText *pt; 			//	http://root.cern.ch/root/html/TPaveText.html
	TImage *resim;			//	http://root.cern.ch/root/html/TImage.html	
	TLegend *leg_hist;		//	http://root.cern.ch/root/html/TLegend.html
	TCanvas *cnvs1,*cnvs2,*cnvs3;			//	http://root.cern.ch/root/html/TCanvas.html
	TString tmpStr, tmpStrName;			//	http://root.cern.ch/root/html/TString.html

	int numOfBins,	numOfBinsSubrange;
	Double_t	x1,x2,y1,y2;

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
	numOfBins=lorentzVectorMvalues->GetNbinsX();


	cout << "\nNumber of events with exactly 2 electrons = " << n2eleevents << endl;
	cout << "\nNumber of bins in 'lorentzVectorMvalues' = " << numOfBins << endl;

	////    TF1 *gfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
	//TF1 *gfit = new TF1("gfit","gaus",minf,maxf); // Create the fit function
	//     lorentzVectorMvalues->Fit("gfit","LL");	//Log-likelihood turu

	/////////////AAAAAAAAAAAAAAAAAAAAAAAAAAAAA/////////////////////////////
	///// 8.7.2012
	double min1=850;
	double max1=1050;	// fit functions will work in range [min1,max1]
	//	3 FIT FUNCTIONS : GAUSSIAN, POLY OF DEG.1, POLY OF DEG.2
	TF1 *fit1=new TF1("fit1","gaus",min1,max1); // Create the fit function
	TF1 *fitpol1=new TF1("fitpol1","pol1",min1,max1); // Create the fit function 
	TF1 *fitpol2=new TF1("fitpol2","pol2",min1,max1); // Create the fit function

	//////////////////////////////////////////////////////////////////////////////
	//	COMPOSITE FIT FUNCTION : GAUSSIAN + POLYNOMIAL OF DEGREE 1
	// Define the parameter array for the total function "fitnoise1total"
	Double_t par1[5];

	// fit functions for "fitnoise1total"
	TF1 *fitnoise1gaus=new TF1("fitnoise1gaus","gaus",min1,max1);
	TF1 *fitnoise1pol=new TF1("fitnoise1pol","pol1",min1,max1);
	TF1 *fitnoise1total=new TF1("fitnoise1total","gaus(0)+pol1(3)",min1,max1);

	cout<<endl<<"*****\t"<<fitnoise1total->GetName()<<"\t*****";

	cout << endl << "\t*****\t" << fitnoise1gaus->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
	lorentzVectorMvalues->Fit("fitnoise1gaus","0R+");

	cout << endl <<"\t*****\t" << fitnoise1pol->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
	lorentzVectorMvalues->Fit("fitnoise1pol","0R+");

	// GET AND SET PARAMETERS OF "fitnoise1total"
	fitnoise1gaus->GetParameters(&par1[0]);
	fitnoise1pol->GetParameters(&par1[3]);
	fitnoise1total->SetParameters(&par1[0]);

	//////////////////////////////////////////////////////////////////////////////
	//	COMPOSITE FIT FUNCTION : GAUSSIAN + POLYNOMIAL OF DEGREE 2
	// Define the parameter array for the total function "fitnoise2total"
	Double_t par2[6];

	// fit functions for "fitnoise2total"
	TF1 *fitnoise2gaus=new TF1("fitnoise2gaus","gaus",min1,max1);
	TF1 *fitnoise2pol=new TF1("fitnoise2pol","pol2",min1,max1);
	TF1 *fitnoise2total=new TF1("fitnoise2total","gaus(0)+pol2(3)",min1,max1);

	cout<<endl<<"*****\t"<<fitnoise2total->GetName()<<"\t*****";

	cout << endl << "\t*****\t" << fitnoise2gaus->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
	lorentzVectorMvalues->Fit("fitnoise2gaus","0R+");

	cout << endl << "\t*****\t" << fitnoise2pol->GetName() << " 'fit' fonksiyonu ile ilgili DATA" <<"\t*****\t" << endl;
	lorentzVectorMvalues->Fit("fitnoise2pol","0R+");

	// GET AND SET PARAMETERS OF "fitnoise1total"
	fitnoise2gaus->GetParameters(&par2[0]);
	fitnoise2pol->GetParameters(&par2[3]);
	fitnoise2total->SetParameters(&par2[0]);

	//////////////////////////////////////////////////////////////////////////////   

	// SET COLOR OF FIT FUNCTIONS
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
	   8	green
	 */

	/*
	 * 	"0" Do not plot the result of the fit. By default the fitted function is drawn unless the option "N"
above is specified
	"R" Use the range specified in the function range
	"+" Add this new fitted function to the list of fitted functions (by default, the previous function is
deleted and only the last one is kept)
	 */

	//////////////////	 CANVAS #1	//////////////////
	cnvs1=new TCanvas();

	// PERFORM FIT OPERATIONS

	// prepare array of fit functions
	int lenOfArr=5;
	TF1 fitFunctionArray[lenOfArr];
	fitFunctionArray[0]=fit1;
	fitFunctionArray[1]=fitpol1;
	fitFunctionArray[2]=fitpol2;
	fitFunctionArray[3]=fitnoise1total;
	fitFunctionArray[4]=fitnoise2total;
	
	// prepare array of legend headers which identify the fit functions
	TString legendHeaderArray[lenOfArr];
	legendHeaderArray[0]="gaus = p0*exp(-0.5*((x-p1)/p2)^2)";
	legendHeaderArray[1]="pol1 = p0 + p1*x";
	legendHeaderArray[2]="pol2 = p0 + p1*x + p2*x^2";
	legendHeaderArray[3]="gaus+pol1";
	legendHeaderArray[4]="gaus+pol2";

//	TF1 *tmpFitFuncLL;

	for(int i=0;i<lenOfArr;i++)
	{
		cout << endl << "\t*****\t" << fit1->GetName() << " 'fit' fonksiyonu ile ilgili DATA, fit : chi2" <<"\t*****\t" << endl;
		fitFunctionArray[i]->SetLineColor(2);
		tmpStr=fitFunctionArray[i]->GetName();
		lorentzVectorMvalues->Fit(tmpStr,"R+");

		// LOG LIKELIHOOD FIT
		TF1 *tmpFitFuncLL=fitFunctionArray[i]->Clone();
		tmpStr+="_LL";
		tmpFitFuncLL->SetName(tmpStr);
		tmpFitFuncLL->SetLineColor(3);
		lorentzVectorMvalues->Fit(tmpStr,"LLR+");

		lorentzVectorMvalues->Draw();
	
		gPad->Update();
		st = (TPaveStats*)lorentzVectorMvalues->FindObject("stats");
//		st->SetOptStat(2211);

		// PUT PAVETEXT ONTO CANVAS, IT WILL BE UNDER PAVESTATS
		x2=st->GetX2NDC(); y2=st->GetY1NDC(); x1=st->GetX1NDC(); y1=0.0;	// coordinates of "PaveText". put the "PaveText" instance right under "PaveStats".
		y1=y2-.1;		// vertical length of "PaveText" is 0.1 NDC.
		pt = new TPaveText(x1,y1,x2,y2,"brNDC");
		tmpStr="# bins = \t";    			tmpStr+=numOfBins;  				  			pt->AddText(tmpStr);
		tmpStr="rebin Parameter = \t";		tmpStr+=rebinParam;			  					pt->AddText(tmpStr);
		tmpStr="% of events = \t";			tmpStr+=Form("%.1f",percentageOfEvents*100);    pt->AddText(tmpStr);
		pt->Draw();
		
		gPad->Update();
		x2=pt->GetX2NDC();	y2=pt->GetY1NDC();	x1=pt->GetX1NDC()-0.1;	y1=y2-.2;
		leg_hist = new TLegend(x1,y1,x2,y2);
		leg_hist->SetHeader(legendHeaderArray[i]);
		TF1 tmpFit=fitFunctionArray[i];
		tmpStr=fitFunctionArray[i]->GetName();
		leg_hist->AddEntry(fitFunctionArray[i],tmpStr,"l");
		
		tmpStr=tmpFitFuncLL->GetName();
		leg_hist->AddEntry(tmpFitFuncLL,tmpStr,"l");
		leg_hist->Draw();
		
		
		// ZOOM IN AREA OF INTEREST, NAMELY [min1,max1]
		//Like for any other ROOT object derived from TObject, one can use the Clone() function. This makes an identical copy of the original histogram including all associated errors and functions
		TH1F *lorentzVectorMvaluesSubRange=lorentzVectorMvalues->Clone();	

		lorentzVectorMvaluesSubRange->SetName("lorentzVectorMvaluesSubRange");
		tmpStr="with range of interest : (";
		tmpStr+=min1;	tmpStr+=",";	
		tmpStr+=max1;	tmpStr+=")";
		lorentzVectorMvaluesSubRange->SetTitle(tmpStr);
		//lorentzVectorMvaluesSubRange->GetXaxis()->SetRange(min1,max1);
		lorentzVectorMvaluesSubRange->GetXaxis()->SetRange(min1/rebinParam,max1/rebinParam);	// if "Rebin()" is called, then new range must be specified accordingly.

		numOfBinsSubrange=lorentzVectorMvaluesSubRange->GetNbinsX();
		
		//////////////////	 CANVAS #3	//////////////////
		cnvs3 = new TCanvas();
		lorentzVectorMvaluesSubRange->Draw();

		//	PAVETEXT FOR CANVAS #3
		pt->Clear();
		tmpStr="# bins = \t";    			tmpStr+=numOfBinsSubrange;			  			pt->AddText(tmpStr);
		tmpStr="rebin Parameter = \t";		tmpStr+=rebinParam;			  					pt->AddText(tmpStr);
		tmpStr="% of events = \t";			tmpStr+=Form("%.1f",percentageOfEvents*100);    pt->AddText(tmpStr);
		pt->Draw();				// draw pavetext
		//	LEGEND FOR CANVAS #3
		leg_hist->Draw();		// "legend" etiketlerini �iz.

		// IMAGE FILE FOR CANVAS #3, NAMELY CANVAS CONTAINING "lorentzVectorMvaluesSubRange"
		tmpStr=fitFunctionArray[i]->GetName();
		tmpStr+="(";		tmpStr+=min1;	tmpStr+=",";	tmpStr+=max1;	tmpStr+=")_";
		tmpStr+=Form("%.1f",percentageOfEvents*100);
		tmpStr+="%_rebinArg=";	tmpStr+=rebinParam;
		tmpStr+=".png";
		resim->FromPad(cnvs3);
		resim->WriteImage(tmpStr);
		
	}

	//cnvs1->SaveAs(resimName);
	//////////////////	 DONE WITH	CANVAS #1	//////////////////





	//////////////////	 DONE WITH	CANVAS #3	//////////////////

//	//////////////////	 .ROOT FILE		//////////////////
//	TFile *dosya=new TFile("total_mass_of_2_el_events.root","UPDATE");
//
//	// prepare the proper directory
//	//	http://root.cern.ch/drupal/content/subdirectories-and-navigation
//	tmpStr=Form("%g",percentageOfEvents*100);	//http://stackoverflow.com/questions/9628645/use-printf-to-format-floats-without-decimal-places-if-only-trailing-0s
//	tmpStr+="%_of_events_is_studied";
//	gDirectory->mkdir(tmpStr);			// gDirectory = current directory 
//	gDirectory->cd(tmpStr);				// change to proper directory
//
//	tmpStr="rebinParam=";	tmpStr+=rebinParam;
//	gDirectory->mkdir(tmpStr);		
//	gDirectory->cd(tmpStr);		
//
//	// write canvases
//	cnvs1->Write(lorentzVectorMvalues->GetName());
//
//	cnvs3->Write(lorentzVectorMvaluesSubRange->GetName());
//
//	// write "lorentzVectorMvalues"
//	tmpStr=Form("%.1f",percentageOfEvents*100);	
//	tmpStr+="%_rebinArg=";
//	tmpStr+=rebinParam;
//	tmpStr+="_lvm";
//	lorentzVectorMvalues->SetName(tmpStr);
//	lorentzVectorMvalues->Write();
//
//	// write "lorentzVectorMvaluesSubRange"
//	tmpStr=Form("%.1f",percentageOfEvents*100);
//	tmpStr+="%_rebinArg=";
//	tmpStr+=rebinParam;
//	tmpStr+="_(";		tmpStr+=min1;	tmpStr+=",";	tmpStr+=max1;	tmpStr+=")";
//	tmpStr+="_lvmSubRange";
//	lorentzVectorMvaluesSubRange->SetName(tmpStr);
//	lorentzVectorMvaluesSubRange->Write();
//
//	dosya->Close();
}


/*
 * What to put inside pavetext
 * 
 *  percentage of events
 *  current	number of bins
 *  rebin parameter
 * 
 * */

/*
	// add the values of parameters to the legend, quited
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

	//leg_hist->AddEntry(fitnoise1total,"gaus+pol1 "+strGetParams,"l");	// given up, makes the legend unnecessarily wide
 */
