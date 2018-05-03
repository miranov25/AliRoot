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


/*
* TODO: add checks for nonexistent projections @Boris
* TODO: add auto canvas->cd(1) after DivideTPad
*/

#include "AliPainter.h"
#include "AliTMinuitToolkit.h"
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
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <map>
#include <bitset>
#include <cstring>


//typedef std::map<int, std::vector<int>() > axisRangesMap;
std::map<TString, TString> AliPainter::genValues;
std::map<TString, TString> AliPainter::fitValues;
std::map<TString, TString> AliPainter::drawValues;
std::vector<TString> AliPainter::rangesVec;
std::map<TString, Double_t> AliPainter::statValues;

/// \brief Method allow to divide pad according to specify properties.
/// \param division      - division string
///                         <nullptr,vertical,horizontal>[div0,div1b, ...]
///                            divi - specify number of pads in row (resp. column)
///                                 - sharing parameter for axis
///                                 - btlrm  (bottom, left, top, right, middle) middle means rl for horizontal and tb for vertical
///                                 - set margin0 in case specified
///                                 - technically attribute can be added to the object name  axis-sharing=""
///                                 - 1lpx30 - means set for pad left margin equal 30pixels
/// \param classID        - adds classID to name of the pad
/// \param style          - adds style to name of the pad in order to using AliDrawStyle
/// \param pad            - input pad to divide
/// \return               - created pad
/*!
*  #### Example use:
*
*  1. Let's load tutorial data:
*       For this you can use
*       \code
*         TFile::SetCacheFileDir(".");
*         TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
*         TTree *tree = (TTree *) finput.Get("hisPtAll");
*         hisArray = new TObjArray();
*         TList *keys = finput->GetListOfKeys();
*         for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
*           TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
*           hisArray->AddLast(o);
*         }
*        \endcode
*  2. Then create the canvas and divide it:
*       \code
*         TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
*         AliPainter::DivideTPad("<horizontal>[1,1,1,1]", "Canvas41", "", canvasQA);
*         canvasQA->cd(1);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*         AliPainter::DrawHistogram("hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)", hisArray);
*       \endcode
*         In this simple case we have four vertical pads:
*         ![<horizontal>\[1b,1t,1,1\]](AliPainter_cxx_example1.png)
*         In case with horizontal position we will have next plot:
*  3. You can sharing chosen axis (see meanings of flags in description.):
*       \code
*         TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
*         AliPainter::DivideTPad("<vertical>[1r,1l,1r,1l]", "Raw,Error", "", canvasQA);
*         canvasQA->cd(1);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*         AliPainter::DrawHistogram("hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)", hisArray);
*        \endcode
*        ![<vertical>\[1r,1l,1r,1l\]](AliPainter_cxx_example2.png)
 * 4. You can specify more than one pads in one raw:
*       \code
*         TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
*         AliPainter::DivideTPad("<horizontal>[1,3,2,1]", "Raw,Error","",canvasQA);
*         canvasQA->cd(1);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*         AliPainter::DrawHistogram("hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(div=1,class=Mass,dOption=E)", hisArray);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*       \endcode
*       ![<horizontal>\[1,3,2,1\]](AliPainter_cxx_example3.png)
*  5. You can specify special flag "m" if you want to join middle pads in column:
*       \code
*         TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
*         AliPainter::DivideTPad("<vertical>[1,3m,2m,1]", "Raw,Error", "", canvasQA);
*         canvasQA->cd(1);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*         AliPainter::DrawHistogram("hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(div=1,class=Mass,dOption=E)", hisArray);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*       \endcode
*       ![<vertical>\[1,3m,2m,1\]](AliPainter_cxx_example4.png)
*  6. All previous case set chosen margin equal 0, but you can set your own value in absolute units,
*     now we support pixel("px" flag) and standard("st" flag) root values (percent from canvas size):
*       \code
*         TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
*         AliPainter::DivideTPad("<horizontal>[1rpx200,1,1,1]", "Raw,Error", "", canvasQA);
*         canvasQA->cd(1);
*         AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
*         AliPainter::DrawHistogram("hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)", hisArray);
*       \endcode
*       ![<horizontal>\[1rpx200,1,1,1\]](AliPainter_cxx_example5.png)
*/
//TODO: add divFlag for inheritance Mother Class (ClassName) should be add to children (className)=>nameOfObject.class(ClassMother,ClassDaughter) @Boris
//TODO: should we change the places of *pad and *classID? @Marian
//TODO: mb remove "<>" from division? @Marian
TPad *AliPainter::DivideTPad(const char *division, const char *classID, const char *style, TPad *pad, Int_t verbose) {

  if (pad == nullptr) {
    TCanvas *canv = new TCanvas("AliPainterCanvas", "AliPainterCanvas", 1200, 800);
    ::Info("AliPainter::DivideTPad", "Input pad is null. We created new one.  pad = new TPad(\"AliPainter\", \"AliPainter\",0,0,1,1);");
    pad = canv;
  }

  Int_t nPads = 0, nRows = 0, maxNCols = 0;
  Double_t xl, yl, xu, yu, a, b, mValue;
  TString position = "";
  TString tempStr = "";
  TString wMargin = ""; // margin zero
  TString units = ""; // set margin
  TString padName = "";

  TObjArray *padRows = TString(division).Tokenize("<>[],");
  position = padRows->At(0)->GetName();
  nRows = padRows->GetEntries() - 1;
  if (verbose == 4) ::Info("AliPainter::DivideTPad", "Number of rows is %d", nRows);

  for (Int_t iRow = 0; iRow < nRows; iRow++)
    if (maxNCols < TString(padRows->At(iRow + 1)->GetName()).Atoi())
      maxNCols = TString(padRows->At(iRow + 1)->GetName()).Atoi();

    for (Int_t iRow = 0; iRow < nRows; iRow++) {

    tempStr = TString(padRows->At(iRow + 1)->GetName());
    Int_t nCols = TString(tempStr(0, 1)).Atoi();
    if (verbose == 4) ::Info("AliPainter::DivideTPad", "Number of columns in %d row is %d", nRows, nCols);
    wMargin = TString(tempStr(1, 1));
    units = TString(tempStr(2, 2));
    mValue = TString(tempStr(4, tempStr.Length())).Atof();

    for (Int_t iCol = 0; iCol < nCols; iCol++) {
      pad->cd();
      xl = iCol / Double_t(nCols);
      yl = (nRows - iRow - 1) / Double_t(nRows);
      xu = (iCol + 1) / Double_t(nCols);
      yu = (nRows - iRow) / Double_t(nRows);

      if (position == "vertical") {
        a = xl;
        b = yl;
        yl = a;       //yl -> xl
        xl = 1 - b;   //xl -> 1 - yl
        a = xu;
        b = yu;
        yu = a;       //yu -> xu
        xu = 1 - b;   //xu -> 1 - yu
      }

      padName = TString::Format("pad[%d]", nPads);
      if (classID != TString("")) padName = TString::Format("%s.class(%s)", padName.Data(), classID);
      if (style != TString("")) padName = TString::Format("%s.style(%s)", padName.Data(), style);

      TPad *newPad = new TPad(padName.Data(), padName.Data(), xl, yl, xu, yu);
      if (verbose == 4) ::Info("AliPainter::DivideTPad", "New pad created: %s", padName.Data());
      if (position == "vertical") {
        newPad->SetTopMargin(nCols * newPad->GetLeftMargin() / maxNCols);
        newPad->SetBottomMargin(nCols * newPad->GetRightMargin() / maxNCols);
      } else {
        newPad->SetLeftMargin(nCols * newPad->GetLeftMargin() / maxNCols);
        newPad->SetRightMargin(nCols * newPad->GetRightMargin() / maxNCols);
      }

      newPad = AliPainter::SetPadMargin(newPad, position, wMargin, units, mValue, iCol, nCols);
      newPad->Draw();
      nPads++;
      newPad->SetNumber(nPads);
    }
  }

  padName = pad->GetName();
  if (classID != TString("")) padName = TString::Format("%s.class(%s)", padName.Data(), classID);
  if (style != TString("")) padName = TString::Format("%s.style(%s)", padName.Data(), style);
  pad->SetName(padName.Data());
  return pad;
}

///
/// \brief Function parses division string from AliPainter::DivideTPad and sets attributes.
/// \param cPad
/// \param position
/// \param units
/// \param value
TPad *AliPainter::SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols) {
  Float_t rlValue = 0.0, btValue = 0.0;
  if (TString(units) == "px") {
    rlValue = mValue / cPad->GetWw();
    btValue = mValue / cPad->GetWh();
  } else if (TString(units) == "st") {
    rlValue = mValue;
    btValue = mValue;
  } else {
    rlValue = 0.0;
    btValue = 0.0;
  }

  if (TString(wMargin) == "r") cPad->SetRightMargin(rlValue);
  if (TString(wMargin) == "l") cPad->SetLeftMargin(rlValue);
  if (TString(wMargin) == "t") cPad->SetTopMargin(btValue);
  if (TString(wMargin) == "b") cPad->SetBottomMargin(btValue);

  if (TString(position) == "vertical") {
    if (nCols > 1 && TString(wMargin) == "m" && iCol == 0) cPad->SetTopMargin(btValue);  //top pad in the m
    if (nCols > 1 && TString(wMargin) == "m" && iCol > 0 && (iCol + 1) < nCols) { // middle pads in the m
      cPad->SetTopMargin(btValue);
      cPad->SetBottomMargin(btValue);
    }
    if (nCols > 1 && TString(wMargin) == "m" && (iCol + 1) == nCols)
      cPad->SetBottomMargin(btValue); //bottom in the m
    if (TString(wMargin) == "m" && nCols == 1) { // only 1 with m
      cPad->SetLeftMargin(rlValue);
      cPad->SetRightMargin(rlValue);
    }
  } else { //the same for horizontal
    if (nCols > 1 && TString(wMargin) == "m" && iCol == 0) cPad->SetRightMargin(rlValue);
    if (nCols > 1 && TString(wMargin) == "m" && iCol > 0 && (iCol + 1) < nCols) {
      cPad->SetRightMargin(0);
      cPad->SetLeftMargin(0);
    }
    if (nCols > 1 && TString(wMargin) == "m" && (iCol + 1) == nCols) cPad->SetLeftMargin(rlValue);
    if (TString(wMargin) == "m" && nCols == 1) {
      cPad->SetTopMargin(btValue);
      cPad->SetBottomMargin(btValue);
    }
  }
  return cPad;
}

///
/// \param graph
/// \param option
void AliPainter::SetMultiGraphTimeAxis(TMultiGraph *graph, TString option) {
  TAxis *axis = nullptr;
  for (Int_t i = 0; i < graph->GetListOfGraphs()->GetEntries(); i++) {
    TGraph *cGraph = (TGraph *) graph->GetListOfGraphs()->At(i);
    if (option.Contains("X")) axis = cGraph->GetXaxis();
    if (option.Contains("Y")) axis = cGraph->GetYaxis();
    if (axis) {
      axis->SetNdivisions(510, kFALSE);
      axis->SetTimeDisplay(1);
      axis->SetTimeFormat("%d/%m");
    }
  }
}

// TODO: check perfomance and change symbols to character constant @Boris
/// \brief Private method for parsing arguments in AliPainter::DrawHistogram
/// \param exprsn - string with arguments
/// \return - vector of arguments
//void AliPainter::FillAllVals(TString exprsn, TString &hisName, TString &projections, std::vector<TString> &fitOptions, std::vector<TString> &rangesStrings, Int_t verbose) {
  //check for match of brackets
//  if (exprsn.CountChar('(') != exprsn.CountChar(')')) {
//    ::Error("AliPainter::DrawHistogram", "check brackets in %s", exprsn.Data());
//    return;
//  }
//  //save each argument from input expression into vector of string
//  std::vector<TString> atts;
//  TString verbStr = "";
//  Int_t match = 0, startIndex = 0, finishIndex = 0;
//  Bool_t isChange = kFALSE;
//  for (Int_t i = 0; i < exprsn.Length(); i++) {
//    if (exprsn(i) == '(' && match == 0) {
//    match++;
//    startIndex = i;
//    isChange = kTRUE;
//    } else if (exprsn(i) == '(' && match > 0) match++;
//    else if (exprsn(i) == ')' && match == 1) {
//      match--;
//      finishIndex = i;
//    } else if (exprsn(i) == ')' && match > 1) match--;
//      if (match == 0 && isChange) atts.push_back(TString(exprsn(startIndex + 1, finishIndex - startIndex - 1)));
//  }

//  if (verbose == 4) {
//    verbStr = "";
//    for (Int_t i = 0; i < (Int_t) atts.size(); i++)
//      verbStr += atts[i] + "|";
//      ::Info("AliPainter::DrawHistogram", "Input expression - %s was transform into array %s", exprsn.Data(),verbStr.Data());
//  }
//
//  //get name of histogram
//  hisName = exprsn(0, exprsn.Index("(", 0));
//  //get string of projections
//  projections = atts[1];
//  //get options for fitting
//  fitOptions = AliPainter::OptionStringParser(atts[2], "()", 6);
//  if (verbose == 4 && atts[2] != "") {
//    verbStr = "";
//    for (Int_t i = 0; i < fitOptions.size(); i++)
//      verbStr += fitOptions[i] + "|";
//      ::Info("AliPainter::DrawHistogram", "Input fit option - %s was transform into array %s", atts[2].Data(), verbStr.Data());
//  }
//
//  //get ranges for hist looping
//  if (atts[0].CountChar(',') > 5) {
//    ::Error("AliPainter::DrawHistogram()", "AliPainter::ArgsParser error. rangesString has more parameters then expected.");
//    return;
//  }
//  if (atts[0] != "") rangesStrings = AliPainter::RangesParser(atts[0]);
//  if (verbose == 4 && atts[0] != "") {
//    verbStr = "";
//    for (Int_t i = 0; i < (Int_t) rangesStrings.size(); i++)
//      verbStr += rangesStrings[i] + "|";
//      ::Info("AliPainter::DrawHistogram", "Input ranges option - %s was transform into array %s", atts[0].Data(), verbStr.Data());
//  }
//  //get drawOption
//  std::vector<TString> drawOptions = AliPainter::OptionStringParser(atts[3], "[]", 6);
//  if (verbose == 4 && atts[3] != TString("")) {
//    AliPainter::RegisterDefaultOptions();
//    verbStr = "";
//    for (Int_t i = 0; i < drawOptions.size(); i++) {
//      if (drawOptions[i] != TString("")) {
//        AliPainter::PandasOptionParser(drawOptions[i]);
//        verbStr += drawOptions[i] + "|";
//      }
//    }
//    ::Info("AliPainter::DrawHistogram", "Input draw option - %s was transform into array %s", atts[3].Data(), verbStr.Data());
//  }
//}

/// \brief Private method for parsing fit options in AliPainter::DrawHistogram
/// \param fitStr - string with fit options
/// \return array of values from inputOptions
//std::vector<TString> AliPainter::ParseOptionString(const char *option, const char d[2], Int_t defSize, Int_t verbose) {
//  TString optStr = TString(option);
//  Int_t arg = 0, startIndex = 0;
//  std::vector<TString> vecOptions;
//  for (Int_t i = 0; i < defSize + 1; i++) vecOptions.push_back(TString(""));
//  if (optStr == TString("")) {
//    ::Error("AliPainter", "AliPainter::OptionStringParser(%s,%s,%d). Options string should not be empty.",
//            optStr.Data(), d, defSize);
//    return vecOptions;
//  }
//  for (Int_t i = 0; i <= optStr.Length(); i++) {
//    if (arg > defSize) break;
//    if (optStr(i) == ',' || i == optStr.Length()) {
//      vecOptions[arg] = TString(optStr(startIndex, i - startIndex));
//      arg++;
//      startIndex = i + 1;
//    } else if (optStr(i) == d[0]) {
//      i = optStr.Index(d[1], i);
//      continue;
//    }
//  }
//  return vecOptions;
//}

void AliPainter::RegisterDefaultOptions() {
  //TODO: perhaps, we should use some aliases? @Marian
  //basic option
  AliPainter::genValues[TString("hisName")] = TString("");
  AliPainter::genValues[TString("projections")] = TString("");
  //fit options
  AliPainter::fitValues[TString("name")] = TString("");
  AliPainter::fitValues[TString("strategy")] = TString("");
  AliPainter::fitValues[TString("option")] = TString("");
  AliPainter::fitValues[TString("ranges")] = TString("");
  AliPainter::fitValues[TString("initPar")] = TString("");
  AliPainter::fitValues[TString("drawOpt")] = TString("");
  //draw options
  AliPainter::drawValues[TString("div")] = TString("");
  AliPainter::drawValues[TString("xlim")] = TString("");
  AliPainter::drawValues[TString("ylim")] = TString("");
  AliPainter::drawValues[TString("zlim")] = TString("");
  AliPainter::drawValues[TString("class")] = TString("");
  AliPainter::drawValues[TString("style")] = TString("");
  AliPainter::drawValues[TString("drawOpt")] = TString("");
  // stat options
//  AliPainter::optionValues[TString("max")]           = TString("");
//  AliPainter::optionValues[TString("min")]           = TString("");
//  AliPainter::optionValues[TString("mean")]          = TString("");
//  AliPainter::optionValues[TString("rms")]           = TString("");
//  AliPainter::optionValues[TString("median")]        = TString("");
}

/// /brief Subsidiary method for RangesParser
/// \param n - number of dimensions
/// \param result - map with strings of axis ranges
/// \param arr - temprorary array
/// \param res - vector of string with values of ranges
void AliPainter::RangesMapToString(Int_t n, axisRangesMap result, std::vector<TString> arr, std::vector<TString> &res) {
  TString tempStr = "";
  if (n > 0) {
    for (UInt_t i = 0; i < result[result.size() - n].size(); i += 2) {
      arr[result.size() - n] = TString::Format("%s,%s", result[result.size() - n][i].Data(),
                                               result[result.size() - n][i + 1].Data());
      RangesMapToString(n - 1, result, arr, res);
    }
  } else {
    for (UInt_t j = 0; j < arr.size(); j++) {
      tempStr += arr[j];
      if (j != arr.size() - 1) tempStr += ",";
    }
    res.push_back(tempStr);
  }
}

/// \brief Private method for filling result map (returned by  AliPainter::RangesParser())
/// \param range - concrete range value or python-like array
/// \param axisNum - num for parsing
/// \param result - output map
void AliPainter::RangesToMap(TString range, Int_t axisNum, axisRangesMap &result) {
  TString num = "";
  Ssiz_t fromStart1 = 0;
  Int_t i = 0;
  Int_t rangeArr[4] = {0, 0, 1, 0}; // iLow, iHigh, step, delta
  std::vector<TString> vRanges;
  if (!range.IsFloat() && range.CountChar(':') > 0) {
    while (range.Tokenize(num, fromStart1, ":")) {
      rangeArr[i] = num.Atoi();
      i++;
    }
    for (Int_t j = rangeArr[0]; j <= rangeArr[1] - rangeArr[3]; j += rangeArr[2]) {
      vRanges.push_back(TString::Itoa(j, 10));
      vRanges.push_back(TString::Itoa(j + rangeArr[3], 10));
    }
    result[axisNum] = vRanges;
  } else if (range.IsFloat()) {
    result[axisNum].push_back(range);
  }
}
/// \brief Method finds histogram in inputArray and draw specified projection according to properties.
/// \param expresion        - histogram draw expression
///                         - syntax
///                           - histogramName(<axisRanges>)(<projection string>)(<fitting string>)(<drawing string>)
///                           - axisRange: @done
///                             - if integer bin range   - TAxis::SetRange(min, max)
///                             - if float   user range  - TAxis::SetRangeUser(min,max)
///                             - if Expression is empty - do not specify anything
///                           - projectionString: (i0,i1) @done
///                             - new projection created THn his = hisInput->Projection(i0,i1....)
///                             - at minimum one dimension should be specified, maximum 3D
///                           - fitting string: (fitterName,fitOption,range,initialParam)
///                             - fitterName - whatever fitter registered in the list of fitters
///                                             (defined in AliTMinuitToolkit , maybe also support
///                                              root fit functions)
///                             - fitOption - see AliTMinuitToolkit fitOptions
///                                            we should put there checks of correctness of fit
///                                            options
///                             - range - {x0min,x0max,x1min,xm1max,...} in case not specified - range is not set
///                           - intitialParam - {p0,p1;p2,...;ep0,ep1,...;minp0,minp1,...; maxp0,maxp1 ...} errors, min and max are optionals
/// TODO: may be we should also use TVectorD or TMatrix instead enumeration? @Boris
///                           - drawing string: (padDiv, lims, className, dOption)
///                             - padDiv - allows to choose wich pad use for drawing:
///                               - divFlag = 0 - use the same pad for drawing; (default value)
///                               - divFlag = 1 - use differents pads for drawing;
///                             - lims - allows to set limits for specified axis: (not specified by default)
///                               - xlim = [xmin, xmax] - set the minima and maxima for x-axes;
///                               - ylim = [ymin, ymax] - set the minima and maxima for y-axes;
///                               - zlim = [zmin, zmax] - set the minima and maxima for z-axes;
///                             - className - adds name of class to each object of drawing. It need for applying css style - [AliDrawStyle:ApplyCssStyle();](link to docs):
///                               - class = [Raw,Error] - in the end of name of object will add ".class(Raw,Error)"
///                               - class = Raw - also as class=[Raw] will add .class(Raw)
///                               - class = [] - in this case nothing will add to the end of the name of object (default)
///                             - dOption - root standard draw options [see docs](https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html#draw-options)
/// \param histogramArray   - array of input objects
///                         - Object to draw - histogramArray->FindObject(histogramName)
///                         - in case histogramArray is nullptr or histogram not found gROOT->FindObject()
///                           will be used
/// \param pad              - input pad if nullptr will create new one;
/// \param metaData         - array with metadata describing histogram
///                         - for example in the trees we optionally keep metadata (array of TNamed ()tag,value) in the array "metaTable"
///                         - in case not specified -"metaTable" object from the histogramArray used
/// \param keepArray        - array for keeping temporary objects
///
/// \return
/*!
*   #### Example usage:
*    Data preparation:
*       \code
*         TFile::SetCacheFileDir(".");
*         TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
*         TTree *tree = (TTree *) finput.Get("hisPtAll");
*         hisArray = new TObjArray();
*         TList *keys = finput->GetListOfKeys();
*         for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
*           TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
*           hisArray->AddLast(o);
*         }
*        \endcode
*   Behaviour by default:
*    \code
*      TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
*      AliPainter::DrawHistogram("hisK0DMassQPtTgl()(0)()()", hisArray);
*    \endcode
*    Using of ranges option:
*    1. Simple values specified
*    \code
*     TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
*     AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80)(0)()()", hisArray);
*    \endcode
*    2. Arrays of values specified
*     2.1 Default behaviour
*       By default draw use the same pad, but you can change it specify draw option.
*       \code
*         TCanvas *canvasQA1 = new TCanvas("canvasQA", "canvasQA", 1200, 800);
*         AliPainter::DivideTPad("<horizontal>[1]", "Canvas1", "",  canvasQA1);
*         canvasQA1->cd(1);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E,ylim=[0,])", hisArray);
*       \endcode
*     2.2 Using AliPainter::DivideTPad();
*       In case if you don't want to use the same pad you should specify it.
*       Also we recommend to you use AliPainter::DivideTPad() for create a few pads on your canvas:
*       here we add className(Raw) and you can see how it looks like via canvasQA->ls(). You can use this class name for applying your own styles in css file.
*       \code
*         TCanvas *canvasQA1 = new TCanvas("canvasQA", "canvasQA", 1200, 800);
*         AliPainter::DivideTPad("<horizontal>[1,1]", "Canvas1", "",  canvasQA1);
*         canvasQA1->cd(1);
*         AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E,ylim=[0,], div=1)", hisArray);
*       \endcode
*     2.3 Limitations to y-axis
*         You can set the ranges to yaxis with using second parameter of drawOption.
*     2.4 Using of standard root draw options
*    3 Fitting
*     3.1 Using of standard root fiting methods.
*     3.2 Using of your own fitters created by AliTMinuitToolkit.
*/
//TODO: refactor the code @Boris
void AliPainter::DrawHistogram(const char *expression, const TObjArray *histogramArray, TPad *pad, TObjArray *metaData,
                               TObjArray *keepArray, Int_t verbose) {
//  TString exprsn = expression;
//  TString hisName = "";
//  TString projections = "";
//  TString uniqName = "";
//  THn *hisN = nullptr;
//  std::vector<TString> fitOptions;
//  std::vector<TString> rangesString;
//
//  //parsing arguments
//  AliPainter::ArgsParser(exprsn, hisName, projections, fitOptions, rangesString);
//  // check for existing of histogram
//  hisN = (THn *) histogramArray->FindObject(hisName);
//  if (hisN == nullptr)
//    ::Info("AliPainter::DrawHistogram", "%s not found", (const char *) hisName);
//  Int_t nDims = hisN->GetNdimensions();
//  //checks  for number of dimensions
//  if (nDims < projections.CountChar(',') + 1) {
//    ::Error("AliPainter::DrawHistogram", "%s has only %d dimensions", (const char *) hisName, nDims);
//    return;
//  } else if (nDims > 3)
//    ::Info("AliPainter::DrawHistogram", "You try to draw 4d histogram.");
//  //check does the fitter exist if not register defaults:
//  //TODO: how I can delete it if I use standart fitter. We should add clearing of fPredefinedFitters into AliTMinuitToolkit::ClearData() @Marian
//  if (AliTMinuitToolkit::GetPredefinedFitter(fitOptions[0].Data()) == nullptr) {
//    AliTMinuitToolkit::RegisterDefaultFitters();
//    if (verbose == 4)
//      ::Info("AliPainter::DrawHistogram",
//             "AliTMinuitToolkit::GetPredefinedFitter(%s) is nullptr. AliTMinuitToolkit::RegisterDefaultFitters();",
//             fitOptions[0].Data());
//  }
//  //pad->SetTitle(expression);
//  Int_t histCnt = 1;
//  Ssiz_t fromStart;
//  Double_t xMin = fitOptions[4].Atof(), xMax = fitOptions[5].Atof();
//  Double_t yMin, yMax;
//  TString key;
//
//  std::vector<TString> rangeVec;
//  TString range = "";
//  TString rangeString = "";
//  TString fitStr = "";
//  TString drawString = "";
//  TH1D *hisArray[histCnt];
//  TVirtualPad *cCanvas;
//  TLegend *legend[histCnt];
//  if (pad == nullptr && gPad == nullptr) {
//    ::Error("AliPainter::DrawHistogram", "TPad object doesn't exist.");
//    return;
//  }
//  if (pad == nullptr) cCanvas = gPad->GetMother();
//  if (pad != nullptr && !pad->InheritsFrom("TCanvas")) cCanvas = pad->GetMother();
//  TPad *nextPad;
//
//  if (!rangesString.empty()) histCnt = (Int_t) rangesString.size();
//  if (verbose == 4) ::Info("AliPainter::DrawHistogram", "Count of the histograms is %d", histCnt);
//  //TODO: refactor it. We already have such strings. Let take it from ArgsParser(): @Boris see TODO in AliPainterTest() ArgsPareserTest()
//  for (Int_t i = 0; i < fitOptions.size() && fitOptions[0] != ""; i++) {
//    fitStr += fitOptions[i];
//    if (i != fitOptions.size() - 1) fitStr += ",";
//  }
//  std::map<TString, TString>::iterator it;
//  for (it = optionValues.begin(); it != optionValues.end(); ++it) {
//    if (it->second != TString("")) drawString = it->first + TString("=") + it->second + TString(",");
//    if (it->second != TString("") && it->first == TString("dOption"))
//      drawString = it->first + TString("=") + it->second;
//  }
//  for (Int_t j = 0; j < histCnt; j++) {
//    if (!rangesString.empty()) {
//      rangeString = rangesString[j].Data();
//      fromStart = 0;
//      rangeVec.clear();
//      while (rangesString[j].Tokenize(range, fromStart, ",")) rangeVec.push_back(range);
//      for (Int_t i = 0; i < rangeVec.size(); i += 2) {
//        if (rangeVec[i].CountChar('.') > 0 || rangeVec[i + 1].CountChar('.') > 0) {
//          hisN->GetAxis(i / 2)->SetRangeUser(rangeVec[i].Atof(), rangeVec[i + 1].Atof());
//          if (verbose == 4)
//            ::Info("AliPainter::DrawHistogram", "his->GetAxis(%d)->SetRangeUser(%s,%s);", i / 2, rangeVec[i].Data(),
//                   rangeVec[i + 1].Data());
//        } else {
//          hisN->GetAxis(i / 2)->SetRange(rangeVec[i].Atoi(), rangeVec[i + 1].Atoi());
//          if (verbose == 4)
//            ::Info("AliPainter::DrawHistogram", "his->GetAxis(%d)->SetRange(%s,%s);", i / 2, rangeVec[i].Data(),
//                   rangeVec[i + 1].Data());
//        }
//      }
//    }
//    //fixme: such names don't work with AliDrawStyle::ApplyCssStyle() @Boris
//    //uniqName = TString::Format("%s(%s)(%s)(%s)(%s)[%d]", hisName.Data(), rangeString.Data(), projections.Data(), fitStr.Data(), drawString.Data(), j).Data();
//    uniqName = TString::Format("%s[%d]", hisName.Data(), j).Data();
//
//    if (optionValues["class"] != TString("") && optionValues["class"].CountChar('[') == 0)
//      uniqName = TString::Format("%s.class(%s)", uniqName.Data(), optionValues["class"].Data()).Data();
//    else if (optionValues["class"] != TString("") && optionValues["class"].CountChar('[') > 0)
//      uniqName = TString::Format("%s.class(%s)", uniqName.Data(), TString(optionValues["class"](1,
//                                                                                                optionValues["class"].Length() -
//                                                                                                2)).Data()).Data();
//    // TH1
//    if (projections.CountChar(',') + 1 == 1) {
//      hisArray[j] = hisN->Projection(projections.Atoi());
//      hisArray[j]->SetName(uniqName);
//      hisArray[j]->SetTitle(uniqName);
//      if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->SetName(%s)", j, uniqName.Data());
//      if (keepArray != nullptr) {
//        keepArray->AddLast((TObject *) hisArray[j]);
//        if (verbose == 4) ::Info("AliPainter::DrawHistogram", "keepArray->AddLast((TObject*) %s)", uniqName.Data());
//      }
//      if (TString(optionValues["ylim"](1, optionValues["ylim"].Index(",") - 1)) != TString("")) {
//        key = TString(optionValues["ylim"](1, optionValues["ylim"].Index(",") - 1));
//        auto kit = AliPainter::statValues.find(key);
//        if (kit == AliPainter::statValues.end()) yMin = key.Atof();
//        else yMin = AliPainter::statValues[key];
//        hisArray[j]->SetMinimum(yMin);
//        if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->SetMinimum(%f)", j, yMin);
//      }
//      if (TString(optionValues["ylim"](optionValues["ylim"].Index(",") + 1,
//                                       optionValues["ylim"].Index("]") -
//                                       optionValues["ylim"].Index(",") - 1)) != TString("")) {
//        key = TString(optionValues["ylim"](optionValues["ylim"].Index(",") + 1,
//                                           optionValues["ylim"].Index("]") -
//                                           optionValues["ylim"].Index(",") - 1));
//        auto kit = AliPainter::statValues.find(key);
//        if (kit == AliPainter::statValues.end()) yMax = key.Atof();
//        else yMax = AliPainter::statValues[key];
//        hisArray[j]->SetMaximum(yMax);
//        if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->SetMaximum(%f)", j, yMax);
//      }
//      if (optionValues["div"] == TString("1")) {
//        if (pad != nullptr) {
//          if (pad->InheritsFrom("TCanvas")) pad->cd(j + 1);
//          else pad->cd();
//        }
//        hisArray[j]->SetStats(kFALSE);
//        hisArray[j]->Draw(optionValues["dOption"].Data());
//        if (pad != nullptr) {
//          if (!pad->InheritsFrom("TCanvas")) {
//            if (cCanvas->GetListOfPrimitives()->After(pad) != nullptr) {
//              nextPad = (TPad *) cCanvas->GetListOfPrimitives()->After(pad);
//              pad = nextPad;
//            }
//          }
//        }
//        if (verbose == 4)
//          ::Info("AliPainter::DrawHistogram", "hisArray[%d]->Draw(%s)", j, optionValues["dOption"].Data());
//        if (pad == nullptr) {
//          if (cCanvas->GetListOfPrimitives()->After(gPad) != nullptr) {
//            nextPad = (TPad *) cCanvas->GetListOfPrimitives()->After(gPad);
//            nextPad->cd();
//          }
//        }
//      } else {
//        //fixme: here I should use more general way for check of empty pad @Boris
//        if (gPad->GetListOfPrimitives()->GetEntries() > 1 ||
//            !gPad->GetListOfPrimitives()->At(0)->InheritsFrom("TFrame")) {
//          hisArray[j]->SetStats(kFALSE);
//          hisArray[j]->SetBit(TH1::kNoTitle);
//          hisArray[j]->Draw(TString::Format("SAME%s", optionValues["dOption"].Data()).Data());
//          if (verbose == 4)
//            ::Info("AliPainter::DrawHistogram", "hisArray[%d]->Draw(\"SAME%s\")", j,
//                   optionValues["dOption"].Data());
//        } else {
//          hisArray[j]->SetStats(kFALSE);
//          hisArray[j]->SetBit(TH1::kNoTitle);
//          hisArray[j]->Draw(optionValues["dOption"].Data());
//          if (verbose == 4)
//            ::Info("AliPainter::DrawHistogram", "hisArray[%d]->Draw(%s)", j, optionValues["dOption"].Data());
//        }
//      }
//      if (fitOptions[0] != TString("")) {
//        //TODO: add check for standard names and change function to gaus and add standart root fitters @Boris
//        if (AliTMinuitToolkit::GetPredefinedFitter(fitOptions[0].Data()) != nullptr) {
//          if (verbose == 4)
//            ::Info("AliPainter::DrawHistogram", "AliTMinuitToolkit::Fit(%s)", fitStr.Data());
//          AliTMinuitToolkit::Fit(hisArray[j], fitOptions[0], fitOptions[1], fitOptions[2], fitOptions[3],
//                                 xMin, xMax);
//        } else {
//          if (verbose == 4)
//            ::Info("AliPainter::DrawHistogram",
//                   "AliTMinuitToolkit::Fit(%s) - Non supported fitter %s. We will try to use standard root fitter.",
//                   fitStr.Data(), fitOptions[0].Data());
//          if (verbose == 4) ::Info("AliPainter::DrawHistogram", " hisArray[%d]->Fit(%s)", j, fitStr.Data());
//          //BUG: fitOptions[n].Data() doesn't work for ubuntu 16.04, but works fine for macOS.
//          hisArray[j]->Fit(fitOptions[0], fitOptions[1], fitOptions[2], xMin, xMax);
//        }
//      }
//    } else if (projections.CountChar(',') + 1 > 3) {
//      if (keepArray != nullptr) keepArray->AddLast((TObject *) hisN->Projection(projections.CountChar(',') + 1));
//    } else return;
//  }
}
//TODO: add ApplyLimitsValue - call GetLimitsValue and apply it. priority: high @Boris
//Double_t AliPainter::GetLimitValue(TString expression, TPad *cPad) {
//  std::map<TString, TString> numNam;
//  //TODO: use regexp @Boris
//  numNam[TString("max")]    = TString("[0]");
//  numNam[TString("min")]    = TString("[1]");
//  numNam[TString("mean")]   = TString("[2]");
//  numNam[TString("rms")]    = TString("[3]");
//  numNam[TString("median")] = TString("[4]");
//
////TODO: combine objects into arrays and use TMath for calculating of parameters @Boris
//}

void AliPainter::DrawHistogram(const char *expression, THn *hisN, TPad *pad, TObjArray *metaData, TObjArray *keepArray, Int_t verbose) {
  if (hisN == nullptr) {
    ::Error("AliPainter::DrawHistogram", "Input histogram is null");
    return;
  }
    if (pad == nullptr && gPad == nullptr) {
    ::Error("AliPainter::DrawHistogram", "TPad object doesn't exist.");
    return;
  }
  AliPainter::FillAll(expression, verbose);
  //TODO: extend to TH2D, TH3D @Boris
  THn *origHisN = hisN;
  AliPainter::SaveToKeepArray((TObject *) hisN, keepArray, verbose);
  TH1D *tempHis =  (TH1D *) AliPainter::SetProjections(hisN, verbose);
  if (tempHis == nullptr) return;
  AliPainter::SaveToKeepArray((TObject *) tempHis, keepArray, verbose);
  AliPainter::SetFitter(tempHis, verbose);
  AliPainter::SaveToKeepArray((TObject *) tempHis, keepArray, verbose);
  AliPainter::SetDrawingOptions(tempHis, verbose);
  AliPainter::SaveToKeepArray((TObject *) tempHis, keepArray, verbose);


  tempHis->Draw();


}
//TODO: think how to combine  AliPainter::ParseString and AliPainter::ParseOptionString
///
/// \param inpString
/// \param sep - separator
/// \param verbose
/// \return - array of input strings
std::vector<TString> AliPainter::ParseString(const char *inpString, const char *sep, Int_t verbose) {
  std::vector<TString> atts;
  TString exprsn = TString(inpString);
  if (exprsn.CountChar('(') != exprsn.CountChar(')')) {
    ::Error("AliPainter::DrawHistogram", "check brackets in %s", exprsn.Data());
    return atts;
  }

  atts.push_back(exprsn(0, exprsn.Index("(", 0)));

  TString verbStr = "";
  Int_t match = 0, startIndex = 0, finishIndex = 0;
  Bool_t isChange = kFALSE;
  for (Int_t i = 0; i < exprsn.Length(); i++) {
    if (exprsn(i) == sep[0] && match == 0) {
    match++;
    startIndex = i;
    isChange = kTRUE;
    } else if (exprsn(i) == sep[0] && match > 0) match++;
    else if (exprsn(i) == sep[1] && match == 1) {
      match--;
      finishIndex = i;
    } else if (exprsn(i) == sep[1] && match > 1) match--;
      if (match == 0 && isChange) atts.push_back(TString(exprsn(startIndex + 1, finishIndex - startIndex - 1)));
  }

  if (verbose == 4) {
    TString infoString = "";
    for (auto const& value: atts)
      infoString += value + "|";
    ::Info("AliPainter::ParseString", "Input string \"%s\" was parsed to %s", inpString, infoString.Data());
  }

  return atts;
}

///
/// \param inputExpr
/// \param defSize - size of output vector
/// \param sep - separator
/// \param ignoreBrackets - inside this bracket separators will be ignored
/// \param verbose
/// \return - array of input strings
std::vector<TString> AliPainter::ParseOptionString(const char *inputExpr, Int_t defSize, const char sep, const char ignoreBrackets[2], Int_t verbose) {
  std::vector<TString> atts;
  if (std::strncmp(inputExpr,"",1) == 0) {
    atts.push_back(TString(""));
    return atts;
  }
  TString inputTStr = TString(inputExpr);
  auto arg = 0, startIndex = 0;
  if (defSize == 0) defSize = inputTStr.CountChar(sep) + 1;
  if (defSize == 0) atts.push_back(inputTStr);
  for (auto i = 0; i < defSize; i++) atts.push_back(TString(""));

  for (auto i = 0; i <= inputTStr.Length(); i++) {
    if (arg >= defSize) {
      ::Warning("AliPainter::ParseString", "AliPainter::ParseString(\"%s\", \'%c\', \"%s\") more than %d parameters in input string", inputExpr, sep, ignoreBrackets, defSize);
      break;
    }
    if (inputTStr(i) == TString(sep) || i == inputTStr.Length()) {
      atts[arg] = TString(inputTStr(startIndex, i - startIndex));
      arg++;
      startIndex = i + 1;
    } else if (inputTStr(i) == TString(ignoreBrackets)[0]) {
      i = inputTStr.Index(ignoreBrackets[1], i);
      continue;
    }
  }

  if (verbose == 4) {
    TString infoString = "";
    for (auto const& value: atts)
      infoString += value + "|";
    ::Info("AliPainter::ParseString", "Input string \"%s\" was parsed to %s", inputExpr, infoString.Data());
  }

  return atts;
}

/// \brief Private method for parsing draw options in AliPainter::DrawHistogram
/// \param drawStr - string with draw options
/// \return array of values from inputOptions
void AliPainter::ParsePandasString(const TString optionsStr, std::map<TString, TString> &optMap, const char sep, const char ignoreBrackets[2], Int_t verbose) {
  std::vector<TString> options = AliPainter::ParseOptionString(optionsStr.Data(), optionsStr.CountChar('=') , sep, ignoreBrackets, verbose);
  for (UInt_t i = 0; i < options.size(); i++) {
    TString optionStr = options[i];
    if (verbose == 4) ::Info("AliPainter::ParsePandasString", "Input string - \"%s\" was parsed to \"%s\"", optionsStr.Data(), optionStr.Data());
    TString key = "";
    TString value = "";
    key = TString(optionStr(0, optionStr.Index("="))).ReplaceAll(" ", "");
    if (verbose == 4) ::Info("AliPainter::ParsePandasString", "From string - \"%s\" key is \"%s\"", optionStr.Data(), key.Data());
    value = TString(optionStr(optionStr.Index("=") + 1, optionStr.Length())).ReplaceAll(" ", "");
    if (verbose == 4) ::Info("AliPainter::ParsePandasString", "From string - \"%s\" value is \"%s\"", optionStr.Data(), value.Data());
    if (optMap.find(key) == optMap.end()) {
      TString defaultKeys = "";
      for (auto it=optMap.begin(); it!=optMap.end(); ++it)
        defaultKeys += it->first + ",";
      ::Error("AliPainter::DrawHistogram",
              "key \"%s\" not found in the list of default keys: %s", key.Data(), defaultKeys.Data());
      return;
    }
    optMap[key] = value;
  }
}

///
/// \param expression
/// \param verbose
void AliPainter::FillAll(const char *expression, Int_t verbose) {

  std::vector<TString> initAtts = AliPainter::ParseString(expression, "()", verbose);

  AliPainter::RegisterDefaultOptions();
  AliPainter::genValues["hisName"] = initAtts[0];
  AliPainter::genValues["projections"] =initAtts[2];
  AliPainter::ParseRanges(initAtts[1], verbose);
  AliPainter::ParsePandasString(initAtts[3], AliPainter::fitValues, ',', "()", verbose);
  AliPainter::ParsePandasString(initAtts[4], AliPainter::drawValues, ',', "[]", verbose);

  if (verbose == 4) {
    TString rangesString = "";
    //TODO: don't use c++11 constructions
    for (auto const& value: AliPainter::rangesVec)
      rangesString += value + ",";
    TString genString = "";
    for (auto it=AliPainter::genValues.begin(); it!=AliPainter::genValues.end(); ++it)
      genString += "\"" + it->first + "\"" + "=" "\"" + it->second + "\"" + ",";
    TString fitString = "";
    for (auto it=AliPainter::fitValues.begin(); it!=AliPainter::fitValues.end(); ++it)
      fitString += "\"" + it->first + "\"" + "=" "\"" + it->second + "\"" + ",";
    TString drawString = "";
    for (auto it=AliPainter::drawValues.begin(); it!=AliPainter::drawValues.end(); ++it)
      drawString += "\"" + it->first + "\"" + "=" "\"" + it->second + "\"" + ",";

    ::Info("AliPainter::FillAll", "Input string \"%s\" was parsed to: \n general             values: \"%s\"; \n ragnes: \"%s\"; \n fit values: \"%s\"; \n draw                values: \"%s\";", expression, genString.Data(), rangesString.Data(),             fitString.Data(), drawString.Data());
  }
}

// TODO - now we have 2 steps for parsing of ranges option to array of string with simple range values. 1. Initial string into map. 2. Map into array fo strings. First of all we can use vector of vector instead map, and the second we should think how we can avoid this 2 intermediate steps with map. @Boris
/// \brief method for parsing ranges options in AliPainter::DrawHistogram
/// \param ranges - input string with values of ranges
/// \param verbose
void AliPainter::ParseRanges(const TString ranges, Int_t verbose) {
  AliPainter::rangesVec.clear();
  if (std::strncmp(ranges.Data(),"",1) == 0)  {
    AliPainter::rangesVec.push_back(TString(""));
    return;
  }
  TString range = "";
  Ssiz_t fromStart0 = 0;
  Int_t axisNum = 0, numCnt = 1;
  axisRangesMap result;
  while (ranges.Tokenize(range, fromStart0, ",")) {
    if (range.IsFloat()) {
      if (numCnt % 2 == 0) axisNum--;
      numCnt++;
    }
    RangesToMap(range, axisNum, result);
    axisNum++;
  }
  std::vector<TString> arr;
  for (Int_t i = 0; i < (Int_t) result.size(); i++) {
    arr.push_back("");
  }
  std::vector<TString> res;
  AliPainter::RangesMapToString((Int_t) result.size(), result, arr, res);
  AliPainter::rangesVec = res;
}

TObject *AliPainter::SetProjections(THn *inHisN, Int_t verbose) {

    std::vector<TString> proj = AliPainter::ParseOptionString(AliPainter::genValues["projections"].Data());
    TObject *res;
    if (proj.size() == 3) res = (TObject *) inHisN->Projection(proj[0].Atoi(), proj[1].Atoi(), proj[2].Atoi());
    else if (proj.size() == 2) res = (TObject *) inHisN->Projection(proj[0].Atoi(), proj[1].Atoi());
    else if (proj.size() == 1) res = (TObject *) inHisN->Projection(proj[0].Atoi());
    else {
      ::Error("AliPainter::SetProjections", "Number of projections not in {1,2,3}");
      res = nullptr;
    }
    return res;
}

template <typename T>
  void AliPainter::SetFitter(T *&inHis, Int_t verbose) {
  if (AliPainter::fitValues["name"] == TString()) {
    if (verbose == 4)
      ::Info("AliPainter::SetFitter", "Fitter was not specified.");
    return;
  }
  std::vector<TString> fitRanges = AliPainter::ParseOptionString(AliPainter::fitValues["ranges"]);
  Double_t xMin = TString(AliPainter::fitValues["ranges"](1,  AliPainter::fitValues["ranges"].Index(",") - 1)).Atof();
  Double_t xMax = TString(AliPainter::fitValues["ranges"](AliPainter::fitValues["ranges"].Index(",")+1,\
                          AliPainter::fitValues["ranges"].Index(")") - AliPainter::fitValues["ranges"].Index(",")-1)).Atof();
  TString fitStr = "";
  if (verbose == 4) {
    for (auto it = AliPainter::fitValues.begin(); it != AliPainter::fitValues.end(); ++it)
      fitStr += "\"" + it->first + "\"" + "=" "\"" + it->second + "\"" + ",";
  }

  if (AliTMinuitToolkit::GetPredefinedFitter(AliPainter::fitValues["name"].Data()) != nullptr) {
    if (verbose == 4) {
      ::Info("AliPainter::DrawHistogram", "AliTMinuitToolkit::Fit(%s)", fitStr.Data());
    }
    AliTMinuitToolkit::Fit(inHis, AliPainter::fitValues["name"], AliPainter::fitValues["strategy"],\
                           AliPainter::fitValues["option"], AliPainter::fitValues["drawOpt"],\
                           xMin, xMax);
  }
  else {
    if (verbose == 4) {
      ::Info("AliPainter::DrawHistogram", \
             "AliTMinuitToolkit::Fit(%s) - Non supported fitter %s. We will try to use standard root fitter.", \
             fitStr.Data(), AliPainter::fitValues["name"].Data());
      ::Info("AliPainter::DrawHistogram", "%s->Fit(%s)", inHis->GetName(), fitStr.Data());
    }
    inHis->Fit(AliPainter::fitValues["name"],  AliPainter::fitValues["option"],  AliPainter::fitValues["drawOpt"], xMin, xMax);
  }

}

void AliPainter::SaveToKeepArray(TObject *obj, TObjArray *&keepArray, Int_t verbose) {
  if (keepArray != nullptr) {
    keepArray->AddLast(obj);
    if (verbose == 4) ::Info("AliPainter::SaveToKeepArray", "Object %s saved to keepArray %s", obj->GetName(), keepArray->GetName());
    return;
  }
  if (verbose == 4) ::Info("AliPainter::SaveToKeepArray", "KeepArray is null");
}

void AliPainter::SaveToKeepArray(TObjArray *objArr, TObjArray *&keepArray, Int_t verbose) {
  if (keepArray != nullptr) {
    for (auto i = 0; objArr->GetEntriesFast(); i++) {
      keepArray->AddLast(objArr->At(i));
      if (verbose == 4)
        ::Info("AliPainter::SaveToKeepArray", "Object %s saved to keepArray %s", objArr->At(i)->GetName(), keepArray->GetName());
    }
    return;
  }
  if (verbose == 4) ::Info("AliPainter::SaveToKeepArray", "KeepArray is null");
}

template <typename T>
  void AliPainter::SetDrawingOptions(T *&inHis, Int_t verbose) {
  if (AliPainter::drawValues[TString("class")] != TString())
    inHis->SetName(TString::Format("%s.class(%s)", inHis->GetName(), AliPainter::drawValues[TString("class")].Data()));
  if (AliPainter::drawValues[TString("style")] != TString())
    inHis->SetName(TString::Format("%s.style(%s)", inHis->GetName(), AliPainter::drawValues[TString("style")].Data()));
  inHis->SetOption(AliPainter::drawValues[TString("drawOpt")].Data());

}
//TODO: change global maps to local
template <typename T>
TObjArray *AliPainter::SetRanges(T *inHis, Int_t verbose) {
  std::vector<TString> rangeVec;

  for (UInt_t i = 0; i < AliPainter::rangesVec.size(); i++) {
    rangeVec.clear();
    rangeVec = AliPainter::ParseOptionString(AliPainter::rangesVec[i].Data());
    for (UInt_t i = 0; i < rangeVec.size(); i += 2) {
        if (rangeVec[i].CountChar('.') > 0 || rangeVec[i + 1].CountChar('.') > 0) {
          hisN->GetAxis(i / 2)->SetRangeUser(rangeVec[i].Atof(), rangeVec[i + 1].Atof());
          if (verbose == 4)
            ::Info("AliPainter::DrawHistogram", "his->GetAxis(%d)->SetRangeUser(%s,%s);", i / 2, rangeVec[i].Data(),
                   rangeVec[i + 1].Data());
        } else {
          hisN->GetAxis(i / 2)->SetRange(rangeVec[i].Atoi(), rangeVec[i + 1].Atoi());
          if (verbose == 4)
            ::Info("AliPainter::DrawHistogram", "his->GetAxis(%d)->SetRange(%s,%s);", i / 2, rangeVec[i].Data(),
                   rangeVec[i + 1].Data());
        }
      }
    }
  }
}
