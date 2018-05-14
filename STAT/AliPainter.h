#ifndef ALIPAINTER_H
#define ALIPAINTER_H

/*
 gSystem->AddIncludePath("-I$AliRoot_SRC/STAT");
.L $AliRoot_SRC/STAT/AliPainter.cxx+

TFile::SetCacheFileDir(".");
TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
TTree *tree = (TTree *) finput->Get("hisPtAll");
hisArray = new TObjArray();
TList *keys = finput->GetListOfKeys();
for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
hisArray->AddLast(o);
}

THn *hisN = hisArray->FindObject("hisK0DMassQPtTgl");
AliPainter::DrawHistogram("hisName()(1)(name=gaus,option=W)(class=Raw,drawOpt=E)", hisN, 0, 0, 0,4)

*/
/// \ingroup STAT
/// \class AliPainter
/*!
* \brief Class for generating QA reports
*  See the documentation in describing of functions.
* \author  <a href="marian.ivanov@cern.ch">Marian  Ivanov</a> , <a href="boris.rumyantsev@cern.ch">Boris Rumyantsev</a>
*/

#include <vector>
#include <map>
#include "TObjArray.h"
#include "TString.h"
#include "TPad.h"
#include "THn.h"
#include "TMultiGraph.h"
#include "TFormula.h"
#include "TObject.h"

class AliPainter : public TObject {
  public:
    static TPad *DivideTPad(const char *division, const char *classID="", const char *style="", TPad *pad=nullptr, Int_t verbose=0);
    static void SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
    static TPad *SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols);
    static TObjArray *PrepareHistogram(const char *expression, THn *hisN, TObjArray *&keepArray, TObjArray *metaData=nullptr, Int_t verbose=0);
    //TODO: hisN should be first argument (in future could be include in THnBase)
    static void DrawHistogram(const char *expression, THn *hisN, TPad *pad=nullptr, TObjArray *keepArray=nullptr, TObjArray *metaData=nullptr, Int_t verbose=0);
    static void DrawHistogram(const char *expression, const TObjArray *histogramArray, TPad *pad=nullptr, TObjArray *keepArray=nullptr, TObjArray *metaData=nullptr, Int_t verbose=0);
    static TPad *GetNextPad(TPad *cPad, TPad *tempPad=nullptr, Int_t verbose=0);

    //static Double_t GetLimitValue(TString expression, TPad *cPad);
    //static void ApplyLimitValue();
    typedef std::map<Int_t, std::vector<TString> > axisRangesMap;
    static std::map<TString, TString> drawValues;
    static std::map<TString, TString> fitValues;
    static std::map<TString, TString> genValues;
    static std::vector<TString> rangesVec;
    static void FillAll(const char *expression, Int_t verbose=0);
    static void ParseRanges(const TString ranges, Int_t verbose=0);
    static void SaveToKeepArray(TObject *obj, TObjArray *&keepArray, Int_t verbose=0);
    static void SaveToKeepArray(TObjArray *objArr, TObjArray *&keepArray, Int_t verbose=0);
    static TObjArray *SetRanges(THn *, Int_t verbose=0);
    static TObject *SetProjections(THn *inHis, Int_t verbose=0);
    static Double_t *GetDataArray(TObjArray *hisArr, Long64_t &commonSize, Int_t verbose=0);
    static void SetLimits(TObjArray *&hisArr, Int_t verbose=0);
    static Double_t GetStatVal(Double_t *valuesArray, Long64_t commonSize, const TString statUnit, Int_t verbose=0);
    template <typename T>
      static void SetFitter(T *&inHis, Int_t verbose=0);
    template <typename T>
      static void SetDrawingOptions(T *&inHis, Int_t verbose=0);
    static std::vector<TString> ParseString(const char *iString,  const char *sep="()", Int_t verbose=0);
    static std::vector<TString> ParseOptionString(const char *optionsStr, Int_t defSize=0, const char sep=',', const char ignoreBrackets[2]="()", Int_t verbose=0);
    static void ParsePandasString(const TString optionsStr, std::map<TString, TString> &optMap, const char sep=',', const char ignoreBrackets[2]="[]", Int_t verbose=0);
    static void RegisterDefaultOptions();
    static void RangesToMap (TString range , Int_t axisNum, axisRangesMap &result);
    static void RangesMapToString(Int_t n, axisRangesMap result, std::vector<TString> arr, std::vector<TString> &res);
    ClassDef(AliPainter,1);
//    static TString GetGenValue(TString name) {return genValues[name];};
//    static void SetGenValue(TString name, TString value) {genValues[name] = value;};
//    static TString GetRangesValue(Int_t index) {return rangesVec[index];};
//    static void SetRangesValue(TString values) {rangesVec.push_back(values);};
//    static void SetRangesValues(std::vector<TString> vec) {rangesVec = vec;};
//    static Int_t GetSizeOfRanges() {return (Int_t) rangesVec.size();};

};
#endif
