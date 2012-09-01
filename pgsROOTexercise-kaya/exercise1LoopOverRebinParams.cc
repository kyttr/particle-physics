#include "exercise1withRebinParam.cc"

#ifdef __CINT__

int exercise1LoopOverRebinParams(int startRebinParam, int endRebinParam, int increment, double percentageParam) {	//
	
	TString tmpStr=Form("%g",percentageParam*100);
	tmpStr+="%_total_mass_of_2_el_events.root";
	TFile dosya(tmpStr,"RECREATE");	//  create a new file, if the file already exists it will be overwritten.
	dosya.Close();
	
	while(startRebinParam<=endRebinParam)
	{
//		.x exercise1withRebinParam.cc(startRebinParam);	// does not run
		exercise1withRebinParam(startRebinParam, percentageParam);
		startRebinParam+=increment;
	}
	return 0;
}

int exercise1LoopOverRebinParams(int startRebinParam, int endRebinParam, int increment) {	//
	
	exercise1LoopOverRebinParams(startRebinParam,endRebinParam,increment,1);
	return 0;
}

int exercise1LoopOverRebinParams(int startRebinParam, int endRebinParam) {	//
	
	exercise1LoopOverRebinParams(startRebinParam,endRebinParam,1);
	return 0;
}


#endif

/*
int exercise1LoopOverRebinParams(int startRebinParam, int endRebinParam, int increment, TObjArray &percentageEventsArr) {	//
	
	TFile dosya("total_mass_of_2_el_events.root","RECREATE");	//  create a new file, if the file already exists it will be overwritten.
	dosya.Close();
	
	double tmpPercentage;
	int j=startRebinParam;
	Long64_t arrLen = percentageEventsArr->GetEntries();
	
	// loop over different values of percentage
	for (Long64_t i=0;i<arrLen;i++)
	{
		tmpPercentage=percentageEventsArr->GetEntry(i);
		// for each percentage value, loop over different values of "rebin" parameters
		while(j<=endRebinParam)
			{
		//		.x exercise1withRebinParam.cc(startRebinParam);	// does not run
				exercise1withRebinParam(startRebinParam);
				j+=increment;
			}	
		j=startRebinParam;
	}
	return 0;
}*/