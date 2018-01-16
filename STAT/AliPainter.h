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
#include "vector"


class AliPainter : public TObject{

public:
  static void     DivideTPad(TPad*pad, const char *division, const char *classID);
  static void     SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
  static TObject* DrawHistogram(char *expresion, const TObjArray* histogramArray=NULL, TObjArray *metaData=NULL, TObjArray * keepArray=NULL);
  static TObject* SetHistogramRange(TObject *hisN, TString range);
  static TPad *SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols);

private:
  ClassDef(AliPainter,1);
  static std::vector<TString> ArgsParser(TString exprsn);
  static std::vector<TString> FitOptParser(TString fitStr);
};

#endif
