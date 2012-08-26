#include "exercise1withRebinParam.cc"

#ifdef __CINT__
int exercise1LoopOverRebinParams(int startRebinParam, int endRebinParam, int increment) {	//
	
	TFile dosya("different_lorentzVectorMvalues_histos.root","RECREATE");	//  create a new file, if the file already exists it will be overwritten.
	dosya.Close();
	
	while(startRebinParam<=endRebinParam)
	{
//		.x exercise1withRebinParam.cc(startRebinParam);	// does not run
		exercise1withRebinParam(startRebinParam);
		startRebinParam+=increment;
	}
	return 0;
}

int exercise1LoopOverRebinParams(int startRebinParam, int endRebinParam) {	//
	exercise1LoopOverRebinParams(startRebinParam,endRebinParam,1);
	return 0;
}
#endif