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


class AliPainter : public TObject{

public:
  static void     DivideTPad(TPad*pad, const char *division);
  static void     SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
  static TObject* DrawHistogram(char *expresion, TObjArray* histogramArray=NULL, TObjArray *metaData=NULL,
                                TObjArray * keepArray=NULL);
  static TObject* SetHistogramRange(TObject *hisN, TString range);
public:

private:
  ClassDef(AliPainter,1);
};

#endif
