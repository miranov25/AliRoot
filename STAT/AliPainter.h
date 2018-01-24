#ifndef ALIPAINTER_H
#define ALIPAINTER_H

/// \ingroup STAT
/// \class AliPainter
/*!
\brief AliPainter
\code
/// TODO add documentation
\endcode

  \author marian  Ivanov marian.ivanov@cern.ch, Boris
*/

class TPad;
class TMultiGraph;
#include "TObject.h"
#include <vector>
#include <map>
#include "TObjArray.h"
#include "TString.h"

class AliPainter : public TObject{

public:
  static void     DivideTPad(TPad *pad, const char *division, const char *classID);
  static void     SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
  static TObject* DrawHistogram(char *expresion, const TObjArray* histogramArray, TPad *pad=NULL, TObjArray *metaData=NULL, TObjArray * keepArray=NULL, Int_t verbose=0);
  static TObject* GetHistogram(TObject *hisN, TString range, std::vector<TString> fitOptions);
  static TPad *SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols);
  static void SetKeepArray(TObjArray *kArray) {keepArray = kArray;}
  static TObjArray *GetKeepArray() {return keepArray;}


private:
  static TObjArray *keepArray;
  typedef std::map<int, std::vector<int> > axisRangesMap;
  ClassDef(AliPainter,1);
  static void ArgsParser(TString exprsn, TString &hisName, TString &projections, std::vector<TString>  &fitOptions, std::vector<TString> &rangesStrings, Int_t verbose = 0);
  static std::vector<TString> FitOptParser(TString inputFitStr);
  static std::vector<TString> RangesParser(TString inputRangesStr);
  static void RangesToMap (TString range , Int_t axisNum, axisRangesMap &result);
  static void RangesMapToString(Int_t n, axisRangesMap result, std::vector<TString> arr, std::vector<TString> &res);
};

#endif
