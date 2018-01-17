
#include <TObjString.h>
#include <TString.h>

void AliTreeTPCMCValTest() {
    
  gSystem->Load("libSTAT");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libTRDcalib");
  gSystem->Load("libT0calib");
//  gSystem->Load("libTOFcalib");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");    

gROOT->LoadMacro("$AliPhysics_SRC/PWGPP/TPC/macros/tpcMCValidation.C+");
  
AliExternalInfo *info = new AliExternalInfo(".", "", 0);
    
TTree* mcTree = info->GetTreeMCPassGuess();

//TObjString mcProd;             //variables for reading RD tree
TString dirName;
char* mcProdName;//[1000];

mcTree->GetBranch("prodName")->SetAddress(&mcProdName);


cout<<mcTree->GetEntries()<<endl;
//for (Int_t i=610; i<612; i++) {       //<mcTree->GetEntries()
  mcProdName=TString("LHC15k1a1").Data();
//    mcTree->GetEntry(i);
//    system(TString::Format("mkdir -p ./%s/",mcProdName).Data());
        system("mkdir -p ./LHC15k1a1/");
        system("cd ./LHC15k1a1/");
    cout<<outputDir<<endl;    
    outputDir="LHC15k1a1";
        cout<<outputDir<<endl; 
    InitTPCMCValidation("LHC15k1a1","passMC","LHC15o", "pass3_lowIR_pidfix",0,0);
//}


}
