/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliPainter.h"
#include <iostream>
#include "TPad.h"
#include "TList.h"
#include "TAxis.h"
#include "TPRegexp.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TError.h"
#include "THn.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

///
/// \brief Method allow to divide pad according to specify properties.
///
///
/// \param pad           - input pad to divide
/// \param division      - division string
///                         <NULL,vertical,horizontal>[div0,div1, ...]
///                            divi - specify number of pads in row (resp. column)
///                                 - sharing parameter for axis
///                                 - btlr  (bottom,lefttop,right)
///                                 - set margin0 in case specified
///                                 - technically attribute can be added to the object name  axis-sharing=""
///                            classID
/*!
  #### Example use:
    \code
    AliExternalInfo info("","",0);
    tree = info.GetTree("QA.TPC","LHC15o","pass1");
    TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
    AliPainter::DivideTPad(canvasQA,"<vertical>[1,3m,2m,1b,1t]");
    canvasQA->cd(1);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIP:run","1","25", "1",1,1,5,0),"ap");
     canvasQA->cd(2);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:run","1","25", "1",1,1,5,0),"ap");
    canvasQA->cd(3);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(4);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(5);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(6);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(7);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(8);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    \endcode
*/
void AliPainter::DivideTPad(TPad*pad, const char *division) {
  // divide pads

  Int_t    nPads   = 0, nRows = 0;
  Double_t xl      = 0.0, yl = 0.0, xu = 0.0, yu = 0.0, a = 0.0, b = 0.0;
  TString  pos     = "";
  TString  tempStr = "";
  TString  s       = "";

  TObjArray *padRows = TString(division).Tokenize("[],");
  pos = padRows->At(0)->GetName();
  nRows = padRows->GetEntries() - 1;

  for (Int_t iRow = 0; iRow < nRows; iRow++) {

    tempStr = TString(padRows->At(iRow + 1)->GetName());
    Int_t nCols = TString(tempStr(0, 1)).Atoi();
    s = TString(tempStr(1, 2));

    for (Int_t iCol = 0; iCol < nCols; iCol++) {
      pad->cd();
      xl = iCol / Double_t(nCols);
      yl = (nRows - iRow - 1) / Double_t(nRows);
      xu = (iCol + 1) / Double_t(nCols);
      yu = (nRows - iRow) / Double_t(nRows);

      if (pos == "vertical") {
        a = xl;
        b = yl;
        yl = a;       //yl -> xl
        xl = 1 - b;   //xl -> 1 - yl
        a = xu;
        b = yu;
        yu = a;       //yu -> xu
        xu = 1 - b;   //xu -> 1 - yu
      }

      TPad *newPad = new TPad(Form("pad%d", nPads), Form("pad%d", nPads), xl, yl, xu, yu);

      if (s == "r") newPad->SetRightMargin(0);
      if (s == "l") newPad->SetLeftMargin(0);
      if (s == "t") newPad->SetTopMargin(0);
      if (s == "b") newPad->SetBottomMargin(0);

      if (pos == "vertical") {
        if (nCols > 1 && s == "m" && iCol == 0) newPad->SetTopMargin(0);  //top pad in the m
        if (nCols > 1 && s == "m" && iCol > 0 && (iCol + 1) < nCols) { // middle pads in the m
          newPad->SetTopMargin(0);
          newPad->SetBottomMargin(0);
        }
        if (nCols > 1 && s == "m" && (iCol + 1) == nCols) newPad->SetBottomMargin(0); //bottom in the m
        if (s == "m" && nCols == 1) { // only 1 with m
          newPad->SetLeftMargin(0);
          newPad->SetRightMargin(0);
        }
      }
      else { //the same for horizontal
        if (nCols > 1 && s == "m" && iCol == 0) newPad->SetRightMargin(0);
        if (nCols > 1 && s == "m" && iCol > 0 && (iCol + 1) < nCols) {
          newPad->SetRightMargin(0);
          newPad->SetLeftMargin(0);
        }
        if (nCols > 1 && s == "m" && (iCol + 1) == nCols) newPad->SetLeftMargin(0);
        if (s == "m" && nCols == 1) {
          newPad->SetTopMargin(0);
          newPad->SetBottomMargin(0);
        }
      }
      newPad->Draw();
      nPads++;
      newPad->SetNumber(nPads);
    }
  }
}


///
/// \param graph
/// \param option
void AliPainter::SetMultiGraphTimeAxis(TMultiGraph *graph, TString option){
  TAxis *axis = NULL;
  for (Int_t i=0; i<graph->GetListOfGraphs()->GetEntries(); i++) {
    TGraph *cGraph=(TGraph *) graph->GetListOfGraphs()->At(i);
    if (option.Contains("X")) axis = cGraph->GetXaxis();
    if (option.Contains("Y")) axis = cGraph->GetYaxis();
    if (axis) {
      axis->SetNdivisions(510, kFALSE);
      axis->SetTimeDisplay(1);
      axis->SetTimeFormat("%d/%m");
    }
  }
}

///
/// \brief Method find histogram in inputArray and draw specified projection according to properties.
/// \param expresion        - histogram draw expression
///                         - syntax
///                           histogramName(<axisRange>)(<projection string>)(operation)(option)
///                             axisRange: @done
///                                if integer bin range   - TAxis::SetRange(min, max)
///                                if float   user range  - TAxis::SetRangeUser(min,max)
///                                if Expression is empty - do not specify anything
///                             projectionString (i0,i1) @done
///                                new projection created THn his = hisInput->Projection(i0,i1....)
///                                at minimum one dimension should be specified, maximum 3D
///                             operation
///
/// \param inputArray       - array of input objects
///                         - Object to draw - histogramArray->FindObject(histogramName)
///                         - in case histogramArray is NULL or histogram not found gROOT->FindObject() will be used
/// \param metaData         - array with metadata describing histogram axis
///                         - for example in the trees we optionally keep metadata (array of TNamed ()tag,value) in the array "metaTable"
///                         - in case not specified -"metaTable" object from the histogramArray used
/// \param keepArray        - array to keep temporary objects
///
/// \return
/*!
   #### Example use:
   \code
   TFile::SetCacheFileDir(".");
   TFile * finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/alice/data/2015/LHC15o/pass1/LHC15o.pass1.B1.Bin0/performanceHisto.root", "CACHEREAD");
   TTree* tree=(TTree*) finput.Get("hisPtAll");
   hisArray = new TObjArray();
   TList *keys = finput->GetListOfKeys();
   for (Int_t iKey = 0; iKey < keys->GetEntries(); iKey++) {
    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
   }
   AliPainter::DrawHistogram("hisK0DMassQPtTgl()(0)()(err)", hisArray, NULL, NULL);
  \endcode
 */

TObject* AliPainter::DrawHistogram(char *expression, TObjArray* inputArray, TObjArray *metaData, TObjArray * keepArray){

  // TString operationType[8]={"f-mean","f-rms","f-ltm","f-ltmsigma","f-gmean","f-grms","f-median","f-gmean"};
  TString exprsn  = expression;
  TString hisName = "";
  TString atts[4] = {"", "", "", ""};
  TPRegexp          attPat("[(].*?[)]");
  TString tStr    = "";
  THn     *his        = NULL;

  if (exprsn.CountChar('(') != exprsn.CountChar(')')) {
    ::Error("AliPainter::DrawHistogram","check brackets in %s", expression);
    return NULL;
  }

  hisName = exprsn(0,exprsn.Index("(",0));
  his     = (THn *) inputArray->FindObject(hisName);;

  if (his == NULL) {
    ::Info("AliPainter::DrawHistogram", "%s not found", (const char*)hisName);
    return NULL;
  }

  int startIndex  = 0;
  int finishIndex = 0;

  for (Int_t i = 0; i < 4; i++) {
    finishIndex = exprsn.Index(")",startIndex);
    tStr = exprsn(attPat,startIndex);
    atts[i] = tStr(1, tStr.Length() - 2);
    startIndex = finishIndex + 1;
  }

  Int_t nDims = his->GetNdimensions();

  if (nDims < atts[1].CountChar(',') + 1) {
    ::Error("AliPainter::DrawHistogram", "%s has only %d dimensions", (const char*)hisName, nDims);
    return NULL; //does ::Error already break the program?
  }

  // TH1
  if (atts[1].CountChar(',') + 1 == 1) {
    TH1 *his1 = NULL;
    his1 = his->Projection(atts[1].Atoi());
    if (SetHistogramRange((TObject *) his1,atts[0]) != NULL) his1 = (TH1 *) SetHistogramRange((TObject *) his1,atts[0]);
    his1->Draw();
    delete his;
    return (TObject *) his1;
  }
  // TH2
  else if (atts[1].CountChar(',') + 1 == 2) {
    TH2 *his2 = NULL;
    his2 = his->Projection(TString(atts[1].Tokenize(",")->At(0)->GetName()).Atoi(),
                           TString(atts[1].Tokenize(",")->At(1)->GetName()).Atoi());
    if (SetHistogramRange((TObject *) his2,atts[0]) != NULL) his2 = (TH2 *) SetHistogramRange((TObject *) his2,atts[0]);
    his2->Draw();
    delete his;
    return (TObject *) his2;
  }
    // TH3
  else if (atts[1].CountChar(',') + 1 == 3) {
    TH3 *his3 = NULL;
    his3 = his->Projection(TString(atts[1].Tokenize(",")->At(0)->GetName()).Atoi(),
                           TString(atts[1].Tokenize(",")->At(1)->GetName()).Atoi(),
                           TString(atts[1].Tokenize(",")->At(2)->GetName()).Atoi());
    if (SetHistogramRange((TObject *) his3,atts[0]) != NULL) his3 = (TH3 *) SetHistogramRange((TObject *) his3,atts[0]);
    his3->Draw();
    delete his;
    return (TObject *) his3;
  }
  else if (atts[1].CountChar(',') + 1 > 3) {
    delete his;
    return (TObject *) his->Projection(atts[1].CountChar(',') + 1);
  }
  else {
    return his;
  }

}
///
/// \param hisN   - Object from TH1, TH2, TH3, THn classes;
/// \param range  - only for Xaxis. e.g. (100, 200) or (100.1, 200.2)
TObject *AliPainter::SetHistogramRange(TObject *hisN, TString range) {
  if (range != "") {
  if (hisN->InheritsFrom("TH1")) {
    TH1 *his1 = (TH1 *) hisN;
    if (range.Index(".") > 0)
      his1->GetXaxis()->SetRangeUser(TString(range.Tokenize(",")->At(0)->GetName()).Atof(),
                                    TString(range.Tokenize(",")->At(1)->GetName()).Atof());
    else
      his1->GetXaxis()->SetRange(TString(range.Tokenize(",")->At(0)->GetName()).Atoi(),
                                TString(range.Tokenize(",")->At(1)->GetName()).Atoi());
    return (TObject *) his1;
  }
  else if (hisN->InheritsFrom("TH2")) {
    TH2 *his2 = (TH2 *) hisN;
    if (range.Index(".") > 0)
    his2->GetXaxis()->SetRangeUser(TString(range.Tokenize(",")->At(0)->GetName()).Atof(),
                                   TString(range.Tokenize(",")->At(1)->GetName()).Atof());
    else
    his2->GetXaxis()->SetRange(TString(range.Tokenize(",")->At(0)->GetName()).Atoi(),
                               TString(range.Tokenize(",")->At(1)->GetName()).Atoi());
    return (TObject *) his2;
  }
  else if (hisN->InheritsFrom("TH3")) {
    TH3 *his3 = (TH3 *) hisN;
    if (range.Index(".") > 0)
    his3->GetXaxis()->SetRangeUser(TString(range.Tokenize(",")->At(0)->GetName()).Atof(),
                                   TString(range.Tokenize(",")->At(1)->GetName()).Atof());
    else
    his3->GetXaxis()->SetRange(TString(range.Tokenize(",")->At(0)->GetName()).Atoi(),
                               TString(range.Tokenize(",")->At(1)->GetName()).Atoi());
    return (TObject *) his3;
  }
  else return NULL;

  }
  else return NULL;
}