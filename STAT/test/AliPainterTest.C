/// \ingroup STAT/test
/// \brief  test methods of AliPainter
/// Example usage
/*
\code
.L $AliRoot_SRC/STAT/test/AliPainterTest.C+
AliPainterTest();
root.exe -b -q  $AliRoot_SRC/STAT/test/AliPainterTest.C+ | tee AliPainterTest.log
\endcode
*/
//TODO: if methods in AliPainter.h will be private we should create new class AliPainterTest inherits from Alipainter. @Boris
#include "AliPainter.h"
#include "TError.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "Rtypes.h"
#include "AliDrawStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TObjArray.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include <cstring>

void AliPainterTest_ParseRanges();
void AliPainterTest_ParseString();
void AliPainterTest_ParseOptionString();
void AliPainterTest_ParsePandasString();

void AliPainterTest_DivideTPad();
//void AliPainterTest_SetMultiGraphTimeAxisTest();
void AliPainterTest_DrawHistogram();
//void AliPainterTest_GetLimitValueTest();
//void AliPainterTest_ApplyLimitValueTest();

void AliPainterTest() {
  AliPainterTest_ParseRanges();
  AliPainterTest_ParseString();
  AliPainterTest_ParseOptionString();
  AliPainterTest_ParsePandasString();
  AliPainterTest_DivideTPad();
//  AliPainterTest_SetMultiGraphTimeAxisTest();
  AliPainterTest_DrawHistogram();
//  AliPainterTest_GetLimitValueTest();
//  AliPainterTest_ApplyLimitValueTest();
}

void AliPainterTest_ParseOptionString() {
  auto result = 0;
  TString input="gaus,W,fitFunction(1,2,3),E,10,200";
  std::vector<TString> optValuesHandle;
  std::vector<TString> optValues;
  optValuesHandle.push_back(TString("gaus"));
  optValuesHandle.push_back(TString("W"));
  optValuesHandle.push_back(TString("fitFunction(1,2,3)"));
  optValuesHandle.push_back(TString("E"));
  optValuesHandle.push_back(TString("10"));
  optValuesHandle.push_back(TString("200"));
  optValues = AliPainter::ParseOptionString(input, 6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- IsOK", input.Data());
  }
  input="gaus,,fitFunction(),,10.234";
  optValuesHandle.clear();
  optValues.clear();
  optValuesHandle.push_back(TString("gaus"));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString("fitFunction()"));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString("10.234"));
  optValuesHandle.push_back(TString(""));
  optValues = AliPainter::ParseOptionString(input, 6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- IsOK", input.Data());
  }
  input="div=0,class=Mass,dOption=E,xlim=[0,10000]";
  optValuesHandle.clear();
  optValues.clear();
  optValuesHandle.push_back(TString("div=0"));
  optValuesHandle.push_back(TString("class=Mass"));
  optValuesHandle.push_back(TString("dOption=E"));
  optValuesHandle.push_back(TString("xlim=[0,10000]"));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString(""));
  optValues = AliPainter::ParseOptionString(input,6, ',', "[]");
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- IsOK", input.Data());
  }
  input="";
  optValuesHandle.clear();
  optValues.clear();
  optValues = AliPainter::ParseOptionString(input);
  if (optValues.size() != 1 && std::strncmp(optValues[0].Data(), "", 1) != 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- IsOK", input.Data());
  }
  input = "1;2;3;4";
  optValues.clear();
  optValues = AliPainter::ParseOptionString(input);
  if (optValues.size() != 1 && std::strncmp(optValues[0].Data(), "1;2;3;4", 7) != 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\", 6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\", 6)- IsOK", input.Data());
  }

  input="1,2,3,4,5";

  optValuesHandle.clear();
  optValues.clear();
  optValuesHandle.push_back(TString("1"));
  optValuesHandle.push_back(TString("2"));
  optValues = AliPainter::ParseOptionString(input,2);
  for (Int_t i = 0; i < 2; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- IsOK", input.Data());
  }

}

void AliPainterTest_ParsePandasString() {
  AliPainter::RegisterDefaultOptions();
  TString input="div=0";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["div"] != TString("0")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="class=Mass";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["class"] != TString("Mass")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="drawOpt=E";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["drawOpt"] != TString("E")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="xlim=[10.123,20.435]";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["xlim"] != TString("[10.123,20.435]")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="zlim=[]";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["zlim"] != TString("[]")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }

  input="ylim=";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["ylim"] != TString("")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }

  input="zlim=[], xlim=[10.123,20.435], ylim=, class=Mass";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["zlim"] != TString("[]") && AliPainter::drawValues["xlim"] != TString("[10.123,20.435]") && \
      AliPainter::drawValues["ylim"] != TString("") && AliPainter::drawValues["class"] != TString("Mass")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }

}

void AliPainterTest_ParseRanges() {
  TString input = "10,20";
  AliPainter::ParseRanges(input);

  if (AliPainter::rangesVec[0] != input) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "10.23,20,54,32.45,43,65.34";
  AliPainter::ParseRanges(input);
  if (AliPainter::rangesVec[0] != input) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "0:20:10:10";
   AliPainter::ParseRanges(input);
  std::vector<TString> outPutHandle;
  outPutHandle.push_back("0,10");
  outPutHandle.push_back("10,20");
  auto result = 0;
  if (outPutHandle.size() != AliPainter::rangesVec.size()) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  for (Int_t i = 0; i < (Int_t) outPutHandle.size(); i++) if (outPutHandle[i] != AliPainter::rangesVec[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "10:20:10,0:40:10:20,30,40";
  AliPainter::ParseRanges(input);
  outPutHandle.clear();
  outPutHandle.push_back("10,10,0,20,30,40");
  outPutHandle.push_back("10,10,10,30,30,40");
  outPutHandle.push_back("10,10,20,40,30,40");
  outPutHandle.push_back("20,20,0,20,30,40");
  outPutHandle.push_back("20,20,10,30,30,40");
  outPutHandle.push_back("20,20,20,40,30,40");
  result = 0;
  if (outPutHandle.size() != AliPainter::rangesVec.size()) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  for (Int_t i = 0; i < (Int_t) outPutHandle.size(); i++) if (outPutHandle[i] != AliPainter::rangesVec[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
}

void AliPainterTest_ParseString() {
  auto result=0;
  TString input="hisK0DMassQPtTgl(20,80,0:80:20:20,0,10)(0)(name=gaus,drawOpt=W)(class=Mass,drawOpt=E)";

  std::vector<TString> args = AliPainter::ParseString(input);
  if (args[0] == TString("hisK0DMassQPtTgl") && args[2] == TString("0") && args[1] == TString("20,80,0:80:20:20,0,10") && \
      args[3] == TString("name=gaus,drawOpt=W") && args[4] == TString("class=Mass,drawOpt=E")) ::Info("AliPainterTest","AliPainter::ParseString(\"%s\")- IsOK", input.Data());
  else {
    ::Error("AliPainterTest","AliPainter::ParseString(\"%s\")- FAILED", input.Data());
    return;
  }
}

void AliPainterTest_DivideTPad() {
  TCanvas *canvasQA = new TCanvas("canvasQATest", "canvasQATest", 1200, 800);
  AliPainter::DivideTPad("<horizontal>[1b,1m,1r300px,1lst0.3]", "", "", canvasQA);
  canvasQA->Print("canvasQADivideTPadTest.xml");
  canvasQA->Print("canvasQADivideTPadTestFixed.xml");

  auto nDiff = gSystem->GetFromPipe("diff canvasQADivideTPadTest.xml $AliRoot_SRC/STAT/test/canvasQADivideTPadTestFixed.xml  | wc -l").Atoi();
  if (nDiff - 6 <= 0) {
    ::Info("AliPainterTest","AliPainter::DivideTPad(\"canvasQATest\",\"<horizontal>[1b,1m,1m,1lst0.3]\",\"test\")- IsOK");
  }else{
    ::Error("AliPainterTest","AliDrawStyle::DivideTPad(\"canvasQATest\",\"<horizontal>[1b,1m,1m,1lst0.3]\",\"test\")- FAILED");
  }
}

void AliPainterTest_DrawHistogram() {
  TFile::SetCacheFileDir(".");
  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
  TTree *tree = (TTree *) finput->Get("hisPtAll");
  TObjArray *hisArray = new TObjArray();
  TList *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
  TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
  AliPainter::DivideTPad("<horizontal>[1b,1t,1,1]", "Canvas41", "", canvasQA);
  canvasQA->cd(1);
  AliPainter::DrawHistogram((char *) "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray,0,0,0,4);
  AliPainter::DrawHistogram((char *) "hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray,0,0,0,4);
  AliPainter::DrawHistogram((char *) "hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray,0,0,0,4);
  AliPainter::DrawHistogram((char *) "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)",
                            hisArray,0,0,0,4);
  canvasQA->Print("canvasQADrawHistogramTest.xml");

  auto nDiff = gSystem->GetFromPipe("diff canvasQADrawHistogramTest.xml $AliRoot_SRC/STAT/test/canvasQADrawHistogramTestFixed.xml  | wc -l").Atoi();
  if (nDiff - 6 <= 0) {
    ::Info("AliPainterTest",
           "AliPainterTest::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)\",hisArray)- IsOK");
  } else {
    ::Error("AliPainterTest",
            "AliPainterTest::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)\",hisArray)- FAILED");
  }
}