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


///  ### AliDrawStyle - Class to access drawing styles
///  * Several drawing styles can be registered and used in the same moment
///    * Styles are identified using strings as identifiers
///      * TStyle
///      * MarkerStyle[]  AliDrawStyle::GetMarkerStyle(const char *style, Int_t index);
///      * MarkerColors[] AliDrawStyle::GetMarkerColor(const char *style, Int_t index);
///      * FillColors[]   AliDrawStyle::GetFillColor(const char *style, Int_t index);
///  * Default styles are created  AliDrawStyle::SetDefaults()
///    * default style is based on the fig template -  https://twiki.cern.ch/twiki/pub/ALICE/ALICERecommendationsResultPresentationText/figTemplate.C
///    * users should be able to register their own styles (e.g in macros)
///  * Usage (work in progress)
///    * performance reports -  with styles as a parameter
///    * QA reports
///    * AliTreePlayer, and TStatToolkit
/// \author marian  Ivanov marian.ivanov@cern.ch
///
///  ## Example usage
/*!
\code
  AliDrawStyle::SetDefaults()
  // Style example
  //
  AliDrawStyle::PrintStyles(0,TPRegexp("."));
  AliDrawStyle::ApplyStyle("figTemplate");
  gPad->UseCurrentStyle();  // force current style for current data
  //
  // Standard ALICE latex symbols
  AliDrawStyle::PrintLatexSymbols(0,TPRegexp("."))
  AliDrawStyle::GetLatexAlice("qpt")
  AliDrawStyle::AddLatexSymbol("dphi", "#Delta#it#phi (unit)")
  AliDrawStyle::GetLatexAlice("dphi")
  //
  // Standard ALICE marker/colors arrays
  AliDrawStyle::GetMarkerStyle("figTemplate",0)
  AliDrawStyle::GetMarkerColor("figTemplate",0)
\endcode
*/


#include "AliDrawStyle.h"
#include "TStyle.h"
#include "TError.h"
#include "TPRegexp.h"
#include "TColor.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TObjArray.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>

using namespace std;
//
TString AliDrawStyle::fDefaultTStyleID;                            ///< ID of the default TStyle
TString AliDrawStyle::fDefaultArrayStyleID;                        ///< ID of the default array styles
std::map<TString, TString>  AliDrawStyle::fLatexAlice;
std::map<TString, TStyle*>  AliDrawStyle::fStyleAlice;
std::map<TString, TObjArray*>  AliDrawStyle::fCssStyleAlice;       //
std::map<TString, std::vector<int> > AliDrawStyle::fMarkerStyles;  // PLEASE LEAVE THE UNAESTHETIC SPACE
std::map<TString, std::vector<int> > AliDrawStyle::fMarkerColors;  // IN ORDER TO MAKE IT WORK WITH THE
std::map<TString, std::vector<float> > AliDrawStyle::fMarkerSize;  // NATIVE SLC6 COMPILER!!!
std::map<TString, std::vector<int> > AliDrawStyle::fFillColors;
std::map<TString, std::vector<float> > AliDrawStyle::fLineWidth;
std::map<TString, std::vector<float> > AliDrawStyle::fLineStyle;
std::map<TString, std::vector<float> > AliDrawStyle::fLineColor;

void AliDrawStyle::SetDefaults(){
  AliDrawStyle::RegisterDefaultLatexSymbols();
  AliDrawStyle::RegisterDefaultStyle();
  AliDrawStyle::RegisterDefaultMarkers();
}

/// set AliDrawStyle::SetDefaultStyles
/// \param styleName  - default style to be used class function in case of empty style selection
/// \param arrayName   - default style to be used class function in case of empty array style selection
void AliDrawStyle::SetDefaultStyles(const char * styleName, const char* arrayName){
  fDefaultTStyleID=styleName;
  fDefaultArrayStyleID=arrayName;
}

TStyle* RegisterDefaultStyleFigTemplate(Bool_t grayScale);

/// Latex symbol section

/// \param symbol  -id name of latex symbol
/// \return latex symbol to be used in ROOT latex
TString AliDrawStyle::GetLatexAlice(const char * symbol){
  return  fLatexAlice[symbol];
}


/// Get integer from string at index
/// \param format    -  array string
/// \param index     -  element index
/// \param separator -  array separator
/// TODO: using TString - to be replaced by faster variant with rough pointers (make it faster if possible)
/// Example usage:
/*!
\code
   AliDrawStyle::GetIntegerAt("1:4:8",1,":"); // return 4
   AliDrawStyle::GetIntegerAt("1;4;8",2,";"); // return  8
\endcode
 */
Int_t AliDrawStyle::GetIntegerAt(const char * format, Int_t index, const char * separator ){
  if (format==NULL) return -1;
  if (index<0) return -1;
  index++;
  TString sFormat(format);
  TString token(format);
  Int_t position=0;
  Int_t counter=0;
  while (counter<index){
    if (sFormat.Tokenize(token,position,separator)) {
      counter++;
    }else{
      break;
    }
  }
  if (counter==index) return token.Atoi();
  return -1;
}

/// Get integer from string
/// \param format    -  array string
/// \param index     -  element index
/// \param separator -  array separator
/// TODO: using TString - to be replaced by faster variant with rough pointers (Boris check)
/// Example usage:
/*!
\code
   AliDrawStyle::GetFloatAt("1.1:4.1:8.1",1,":"); // return 4.1
   AliDrawStyle::GetFloatAt("1.1;4.1;8.1",2,";"); // return  8.1
\endcode
*/
Float_t AliDrawStyle::GetFloatAt(const char * format, Int_t index, const char * separator ){
  if (format==NULL) return -1;
  if (index<0) return -1;
  index++;
  TString sFormat(format);
  TString token(format);
  Int_t position=0;
  Int_t counter=0;
  while (counter<index){
    if (sFormat.Tokenize(token,position,separator)) {
      counter++;
    }else{
      break;
    }
  }
  if (counter==index) return token.Atof();
  return -1;
}


// GetMarkerStyle associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker style for given styleName, index
Int_t AliDrawStyle::GetMarkerStyle(const char *style, Int_t index){
  if (AliDrawStyle::fMarkerStyles[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fMarkerStyles[style][index];
}

// GetLineStyle associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker style for given styleName, index
Int_t AliDrawStyle::GetLineStyle(const char *style, Int_t index){
  if (AliDrawStyle::fLineStyle[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fLineStyle[style][index];
}

// GetLineColor associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker style for given styleName, index
Int_t AliDrawStyle::GetLineColor(const char *style, Int_t index){
  if (AliDrawStyle::fLineColor[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fLineColor[style][index];
}

/// GetMarkerColor associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker color for given styleName, index
Int_t AliDrawStyle::GetMarkerColor(const char *style, Int_t index){
  if (AliDrawStyle::fMarkerColors[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fMarkerColors[style][index];
}
/// GetMarkerSize associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker color for given styleName, index
Float_t AliDrawStyle::GetMarkerSize(const char *style, Int_t index){
  if (AliDrawStyle::fMarkerSize[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fMarkerSize[style][index];
}
/// GetFillColor associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return fill color for given styleName, index
Int_t AliDrawStyle::GetFillColor(const char *style, Int_t index){
  if (AliDrawStyle::fFillColors[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fFillColors[style][index];
}
/// GetLineWidth associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return fill color for given styleName, index
Float_t AliDrawStyle::GetLineWidth(const char *style, Int_t index){
  if (AliDrawStyle::fLineWidth[style].size() <= index) {
    return GetFloatAt(style,index);
  }
  return  AliDrawStyle::fLineWidth[style][index];
}



void AliDrawStyle::PrintLatexSymbols(Option_t */*option*/, TPRegexp& regExp){
  //print latex symbols
  typedef std::map<TString,TString>::const_iterator it_type;
  for(it_type iterator = fLatexAlice.begin(); iterator != fLatexAlice.end(); ++iterator) {
    if (regExp.Match(iterator->first.Data())==0) continue;
    std::cout<<iterator->first << " " << iterator->second << "\n";
  }
}

void AliDrawStyle::PrintStyles(Option_t *option, TPRegexp& regExp){
  //print latex symbols
  typedef std::map<TString,TStyle*>::const_iterator it_type;
  for(it_type iterator = fStyleAlice.begin(); iterator != fStyleAlice.end(); ++iterator) {
    if (regExp.Match(iterator->first.Data())==0) continue;
    if (option==NULL) std::cout << iterator->first << " " << iterator->second << "\n";
    if (option!=NULL) {
      iterator->second->Print(option);
      if (TString(option).Contains("dump")) iterator->second->Dump();
    }
  }
}


void AliDrawStyle::ApplyStyle(const char* styleName){
  //
  // Apply registered style
  //
  TStyle * style= fStyleAlice[styleName];
  if (style==NULL){
    ::Error("AliDrawStyle::ApplyStyle","Invalid style %s",styleName);
  }else{
    ::Info("AliDrawStyle::ApplyStyle","%s",styleName);
  }
  if (style) style->cd();
}


void  AliDrawStyle::AddLatexSymbol(const char * symbolName, const char * symbolTitle){
  fLatexAlice[symbolName]=symbolTitle;
}
void  AliDrawStyle::RegisterDefaultLatexSymbols(){
  //
  // Set default AliRoot/Latex/root shortcuts
  //
  fLatexAlice["qpt"]="#it{q}/#it{p}_{T} (GeV/#it{c})^{-1}";
  fLatexAlice["qpt0"]="#it{q}/#it{p}_{T}";
  fLatexAlice["pt"]="#it{p}_{T}  (GeV/#it{c}) ";
  fLatexAlice["pt0"]="#it{p}_{T} ";
  fLatexAlice["sqptmev"]="#sigma_{#it{q}/#it{p}_{T}} (MeV/#it{c})^{-1}";
  fLatexAlice["pbpb502"]="Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} =5.02 TeV";
  fLatexAlice["pp13"]="pp #sqrt{#it{s}} = 13 TeV ";
  fLatexAlice["drphi"]="#Delta_{#it{r#phi}} (cm)";
  fLatexAlice["srphi"]="#sigma_{#it{r#phi}} (cm)";
}

void   AliDrawStyle::RegisterDefaultStyle(){
  //
  fStyleAlice["figTemplate"]=RegisterDefaultStyleFigTemplate(kFALSE);
  fStyleAlice["figTemplateGrey"]=RegisterDefaultStyleFigTemplate(kFALSE);
  //
  TStyle *style=RegisterDefaultStyleFigTemplate(kFALSE);
  style->SetName("figTemplate2");
  style->SetTitleXSize(TMath::Power(2,0.5)*style->GetTitleXSize());
  style->SetTitleYSize(TMath::Power(2,0.5)*style->GetTitleYSize());
  style->SetLabelSize(TMath::Power(2,0.5)*style->GetLabelSize("X"),"X");
  style->SetLabelSize(TMath::Power(2,0.5)*style->GetLabelSize("Y"),"Y");
  style->SetLabelSize(TMath::Power(2,0.5)*style->GetLabelSize("Z"),"Z");
  fStyleAlice["figTemplate2"]=style;
  //
  style=RegisterDefaultStyleFigTemplate(kFALSE);
  style->SetName("figTemplate3");
  style->SetTitleXSize(TMath::Power(3,0.5)*style->GetTitleXSize());
  style->SetTitleYSize(TMath::Power(3,0.5)*style->GetTitleYSize());
  style->SetLabelSize(TMath::Power(3,0.5)*style->GetLabelSize("X"),"X");
  style->SetLabelSize(TMath::Power(3,0.5)*style->GetLabelSize("Y"),"Y");
  style->SetLabelSize(TMath::Power(3,0.5)*style->GetLabelSize("Z"),"Z");
  fStyleAlice["figTemplate3"]=style;



}

void  AliDrawStyle::RegisterDefaultMarkers(){
  //
  // Style source:
  // https://twiki.cern.ch/twiki/pub/ALICE/ALICERecommendationsResultPresentationText/figTemplate.C
  const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7, kBlack, kRed+1 }; // for systematic bands
  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2,kGray+1,  kRed-10 };
  const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
  //
  (fMarkerStyles["figTemplate"])=std::vector<int>(10);
  (fMarkerColors["figTemplate"])=std::vector<int>(10);
  (fMarkerSize["figTemplate"])=std::vector<float>(10);
  (fFillColors["figTemplate"])=std::vector<int>(10);
  (fLineWidth["figTemplate"])=std::vector<float>(10);
  (fLineColor["figTemplate"])=std::vector<float>(10);
  (fLineStyle["figTemplate"])=std::vector<float>(10);
  for (Int_t i=0;i<10; i++){
    (fMarkerStyles["figTemplate"])[i]=markers[i];
    (fMarkerColors["figTemplate"])[i]=colors[i];
    (fMarkerSize["figTemplate"])[i]=1;
    (fFillColors["figTemplate"])[i]=fillColors[i];
    (fLineWidth["figTemplate"])[i]=0.5;
    (fLineStyle["figTemplate"])[i]=i+1;
    (fLineColor["figTemplate"])[i]=colors[i];
  }
  // style inspired by TRD performance paper
  Int_t colorsTRD[12]={0};
  const Int_t markersTRD[]    = {kOpenCircle,kFullCircle, kOpenSquare,kFullSquare, kOpenStar,kFullStar, kOpenDiamond,kFullDiamond, kOpenCross,kFullCross };
  const Float_t markerTRDSize[]    = {1,1, 0.9,0.9, 1.4,1.4, 1.1,1.1, 1.2,1.2 };
  colorsTRD[0]=TColor::GetColor("#0000DD");
  colorsTRD[1]=TColor::GetColor("#00EE00");
  colorsTRD[2]=TColor::GetColor("#FF0000");
  colorsTRD[3]=TColor::GetColor("#00EEDD");
  colorsTRD[4]=TColor::GetColor("#FFEE00");
  colorsTRD[5]=TColor::GetColor("#FF00DD");
  colorsTRD[6]=TColor::GetColor("#9999DD");
  colorsTRD[7]=TColor::GetColor("#99EE99");
  colorsTRD[8]=TColor::GetColor("#FF9999");
  colorsTRD[9]=TColor::GetColor("#66AADD");
  colorsTRD[10]=TColor::GetColor("#AAEE66");
  colorsTRD[11]=TColor::GetColor("#FF66AA");
  (fMarkerStyles["figTemplateTRD"])=std::vector<int>(10);
  (fMarkerColors["figTemplateTRD"])=std::vector<int>(10);
  (fMarkerSize["figTemplateTRD"])=std::vector<float>(10);
  (fFillColors["figTemplateTRD"])=std::vector<int>(10);
  (fLineWidth["figTemplateTRD"])=std::vector<float>(10);
  (fLineStyle["figTemplateTRD"])=std::vector<float>(10);
  (fMarkerStyles["figTemplateTRDPair"])=std::vector<int>(10);
  (fMarkerColors["figTemplateTRDPair"])=std::vector<int>(10);
  (fMarkerSize["figTemplateTRDPair"])=std::vector<float>(10);
  (fFillColors["figTemplateTRDPair"])=std::vector<int>(10);
  (fLineWidth["figTemplateTRDPair"])=std::vector<float>(10);
  (fLineStyle["figTemplateTRDPair"])=std::vector<float>(10);
  for (Int_t i=0; i<10; i++){
    (fMarkerStyles["figTemplateTRD"])[i]=markersTRD[i];
    (fMarkerColors["figTemplateTRD"])[i]=TColor::GetColorDark(colorsTRD[i]);
    (fMarkerSize["figTemplateTRD"])[i]=markerTRDSize[i];
    (fFillColors["figTemplateTRD"])[i]=fillColors[i];
    (fLineWidth["figTemplateTRD"])[i]=0.5;
    //
    (fMarkerStyles["figTemplateTRDPair"])[i]=markersTRD[i];
    (fMarkerColors["figTemplateTRDPair"])[i]=TColor::GetColorDark(colorsTRD[i/2]);
    (fMarkerSize["figTemplateTRDPair"])[i]=markerTRDSize[i];
    (fFillColors["figTemplateTRDPair"])[i]=fillColors[i/2];
    (fLineWidth["figTemplateTRDPair"])[i]=0.5;
  }

}

TStyle*  RegisterDefaultStyleFigTemplate(Bool_t grayPalette) {
  // Style source:
  // https://twiki.cern.ch/twiki/pub/ALICE/ALICERecommendationsResultPresentationText/figTemplate.C
  //
  TStyle * figStyle = new TStyle;
  figStyle->Reset("Plain");
  figStyle->SetOptTitle(0);
  figStyle->SetOptStat(0);

  if(grayPalette) figStyle->SetPalette(8,0);
  else figStyle->SetPalette(1);

  figStyle->SetCanvasColor(10);
  figStyle->SetCanvasBorderMode(0);
  figStyle->SetFrameLineWidth(1);
  figStyle->SetFrameFillColor(kWhite);
  figStyle->SetPadColor(10);
  figStyle->SetPadTickX(1);
  figStyle->SetPadTickY(1);
  figStyle->SetPadBottomMargin(0.15);
  figStyle->SetPadLeftMargin(0.15);
  figStyle->SetHistLineWidth(1);
  figStyle->SetHistLineColor(kRed);
  figStyle->SetFuncWidth(2);
  figStyle->SetFuncColor(kGreen);
  figStyle->SetLineWidth(2);
  figStyle->SetLabelSize(0.045,"xyz");
  figStyle->SetLabelOffset(0.01,"y");
  figStyle->SetLabelOffset(0.01,"x");
  figStyle->SetLabelColor(kBlack,"xyz");
  figStyle->SetTitleSize(0.05,"xyz");
  figStyle->SetTitleOffset(1.25,"y");
  figStyle->SetTitleOffset(1.2,"x");
  figStyle->SetTitleFillColor(kWhite);
  figStyle->SetTextSizePixels(26);
  figStyle->SetTextFont(42);
  figStyle->SetLegendBorderSize(0);
  figStyle->SetLegendFillColor(kWhite);
  figStyle->SetLegendFont(42);
  return figStyle;
}

/// CSS style parsing functions

///
/// \param input        - input string
/// \param propertyName - name of property to find
/// \return             - property propertyName from the input string using CSS like parsing, empty string in case not found
///
TString  AliDrawStyle::GetPropertyValue(TString input, TString propertyName){
  Int_t index0 = input.Index(propertyName.Data());
  if (index0<0) return "";
  Int_t index1 = input.Index(':',index0)+1;
  Int_t index2= input.Index(';',index1)-1;
  if (index2-index1-1<0) return "";
  TString result(input(index1,index2-index1+1));
  return result;
}
///
/// \param input    - input string (CSS record - find proper name in the w3c)
/// \param propertyName  -  name of property to find
/// \param index    - index of value to find
/// \return         - value as a integer  -1 if does not exist, resp , separated value at index index
/*!
 ####  Example use:
  TString input="{\nmarker_style:25,21,22,23; \nmarker_color:1,2,4,5; \n}";
  AliDrawStyle::GetNamedIntegerAt(input,"marker_style",0);  //  return  25
  AliDrawStyle::GetNamedIntegerAt(input,"marker_style",3);  //  return  23
  AliDrawStyle::GetNamedIntegerAt(input,"marker_color",2);  //  return  4
 */
Int_t    AliDrawStyle::GetNamedIntegerAt(TString input, TString propertyName, Int_t index){
  TString  value = AliDrawStyle::GetPropertyValue(input,propertyName);
  Int_t indexStart=0;
  Int_t indexFinish=value.Index(',',indexStart);
  for (Int_t j=0; j<index; j++){
    indexStart=value.Index(',',indexStart)+1;
    indexFinish=value.Index(',',indexStart);
    if (indexStart<0 || index>value.CountChar(',')) return -1;
    if (indexFinish<0) indexFinish = value.Length();
  }
  TString valueAt(value(indexStart, indexFinish - indexStart));
  if (valueAt.IsFloat()) return valueAt.Atoi();
  else return -1;
}

/// \param input    - input string (CSS record - find proper name in the w3c)
/// \param propertyName  -  name of tag to find
/// \param index    - index of value to find
/// \return         - value as a float  -1 if does not exist, resp , separated value at index index
/*!
 ####  Example use:
  TString input="{\nmarker_style:25,21,22,23; \nmarker_color:1,2,4,5; \n}";
  AliDrawStyle::GetNamedIntegerAt(input,"marker_style",0);  //  return  25
  AliDrawStyle::GetNamedIntegerAt(input,"marker_style",3);  //  return  23
  AliDrawStyle::GetNamedIntegerAt(input,"marker_color",2);  //  return  4
 */
Float_t  AliDrawStyle::GetNamedFloatAt(TString input, TString propertyName, Int_t index){
  TString  value = AliDrawStyle::GetPropertyValue(input,propertyName);
  Int_t indexStart=0;
  Int_t indexFinish=value.Index(',',indexStart);
  for (Int_t j=0; j<index; j++){
    indexStart=value.Index(',',indexStart)+1;
    indexFinish=value.Index(',',indexStart);
    if (indexStart<0 || index>value.CountChar(',')) return -1;
    if (indexFinish<0) indexFinish = value.Length();
  }
  TString valueAt(value(indexStart, indexFinish - indexStart));
  if (valueAt.IsFloat()) return valueAt.Atof();
  else return -1;
}

/// Read CSS html like files  (*see also AliRoot modification in CSS)
/// TODO:
/// * proper exception  handling (Boris)
///   * code should not fail
///   * return 0 pointer if inconsistent content
///   * Use ::Error verbosity according debug level
/// * proper CSS comments handling (Boris) @done
/// * include CSS files  (should be included as )
/// \param inputName     - input file to read
/// \param verbose       - specify verbose level for ::error and ::info (Int_t should be interpreted as an bit-mask)
/// \return              - TObjArray  with the pairs TNamed of the CSS <Selector, declaration> or  TObjArray (recursive structure like includes)
TObjArray * AliDrawStyle::ReadCSSFile(const char *  inputName, Int_t verbose){
  //check file exist
  TString inputCSS = gSystem->GetFromPipe(TString::Format("cat %s",inputName).Data());     // I expect this variable is defined
  //remove comments:
  while (inputCSS.Index("*/") > 0){
    inputCSS = inputCSS(0, inputCSS.Index("/*")) + inputCSS(inputCSS.Index("*/")+2, inputCSS.Length());
  }
    //inputCSS.ReplaceAll("\n", ""); we can add this, in the other case ClassName will be like "\n .TH*", but I suppose it doesn't matter.
  TObjArray *tokenArray = inputCSS.Tokenize("{}");   //assuming we can not use {} symbols in the style IDS
  Int_t entries = tokenArray->GetEntries();
  TObjArray *cssArray = new TObjArray(entries / 2);
  for (Int_t i = 0; i < entries; i += 2) {
    if (i + 1 >= entries) continue;
    TString selector = tokenArray->At(i)->GetName();
    TString declaration = tokenArray->At(i + 1)->GetName();
    cssArray->AddLast(new TNamed(selector.Data(), declaration.Data()));
  }
  return cssArray;
}

/// Write cssArray to the file as a plain array (recursive function)
/// \param cssArray    - input css array to write
/// \param outputName  - output file
/// \param cssOut     - output stream ( )
void    AliDrawStyle::WriteCSSFile(TObjArray * cssArray, const char *  outputName, fstream *pCssOut) {
  if (pCssOut == NULL) {
    pCssOut=new fstream;
    pCssOut->open("test.css", ios_base::out|ios_base::trunc);
  }
  fstream &cssOut = *pCssOut;
  for (Int_t i=0;i<cssArray->GetEntries();i++) {
    TObject *object = cssArray->At(i);
    if (object->InheritsFrom("TObjArray")){
      AliDrawStyle::WriteCSSFile((TObjArray*)object, 0, pCssOut);
    }else {
      cssOut << object->GetName();
      cssOut << "{";
      cssOut << object->GetTitle();
      cssOut << "}";
    }
  }
  cssOut<<std::flush;
  if (outputName!=NULL) {
    pCssOut->close();
    delete pCssOut;
  }
}
/// Function to check  match between "CSS" selector and pair of className, objectName
/// \param selectors    - selector ID
/// \param elementName - name of element
/// \param className   - name of class
/// \param objectName  - object name
/// \return            - kTRUE if selector match class name and object name
/// TODO
///   - pre-calculate the tree ? yes, in best practice we should use trees for css parsing. I think later we can implement simple variant.
Bool_t  AliDrawStyle::IsSelected(TString selectors, TString elementName, TString className, TString objectName){
  //TString selectors = "TH1.Status#obj1, TH1.Warning#obj1, TH1.Warning#obj3 \tTGraph#obj1, TGraph.Status#TPC.QA.dcar_posA_1 \tTGraph.Warning#TPC.QA.dcar_posA_2 \tTF1.Status, .Status#obj1, #obj3"

  Bool_t elementCatched;
  Bool_t classCatched;
  Bool_t objectCatched;
  Ssiz_t fromStart=0;
  Ssiz_t fromStart1;

  TString subSelectors;
  TString selector;

  TPRegexp elemPat(".*?[.]");
  TPRegexp classPat("[.].*[#]");
  TPRegexp objPat("[#].*");

  while(selectors.Tokenize(subSelectors,fromStart," \t")){
    fromStart1 = 0;
   while(subSelectors.Tokenize(selector,fromStart1,", ")){
     selector = selector.Strip(TString::kBoth, ' ');
     elementCatched = false;
     classCatched = false;
     objectCatched = false;

     if(selector(objPat) == ""){
        objectCatched = true;
        selector.Append("#anyObjects");
      }
      else if (selector(objPat) == ("#" + objectName)) objectCatched = true;
      if(selector(classPat) == ""){
         classCatched = true;
         selector.Insert(selector.Index("#"), ".anyClasses");
       }
       else if (selector(classPat) == ("." + className + "#")) classCatched = true;
     if(selector(elemPat) == "." || selector(elemPat) == "" || selector(elemPat) == (elementName + ".")) elementCatched = true;
     if (elementCatched && classCatched && objectCatched)  return true;
   }
  }
  return false;
}

///
/// \param styleName
/// \param elementID
/// \param classID
/// \param objectID
/// \return
/*!
\code
AliDrawStyle::SetCssStyle("alirootTestStyle.css",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css",0));
AliDrawStyle::GetProperty("alirootTestStyle.css","marker_size", "TGraph", "", "obj1");                     // "   1,1,1,1" as it should be
AliDrawStyle::GetProperty("alirootTestStyle.css","marker_size", "TGraph", "", "TPC.QA.dcar_posA_1");       // empty as it should
AliDrawStyle::GetProperty("alirootTestStyle.css","marker_size", "TGraph", "Status", "TPC.QA.dcar_posA_1"); // should be "   1,1,1,1"


\endcode
*/
TString AliDrawStyle::GetProperty(const char *styleName, TString propertyName, TString elementID, TString classID, TString objectID){

  if(fCssStyleAlice[styleName] == NULL) return "";

  Int_t entries = fCssStyleAlice[styleName]->GetEntriesFast();
  TString declaration="";
  for(Int_t i = 0; i < entries; i++){
    if(IsSelected(TString(fCssStyleAlice[styleName]->At(i)->GetName()), elementID, classID, objectID)){
      TString  value = GetPropertyValue(fCssStyleAlice[styleName]->At(i)->GetTitle(),propertyName);
      if (value.Length()>0) return value.Strip(TString::kBoth, ' ');
    }
  }
  return "";
}
