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

///   ### AliDrawStyle - Class to access drawing styles
///   Our goal is provide to users comfortable interface for changing style of drawing.
///   We think that best tools for such tasks it's CSS. So it's attempt to create analog of CSS for aliroot.
///   All you need to know about CSS for usage AiDrawStyle is:
///     1. First of all you should to create and fill file with ".css" extension with description of your own styles according to your objects declaration.
///        (See example)[https://github.com/alisw/AliRoot/blob/master/STAT/test/alirootTestStyle.css]
///     2. The syntax of CSS is
///         \image html CSSsyntax.png
///         * selector - define the element (class in root), class (user predefined class (like tag)) via name of object.
///           E.g. TH1 *his = TH1("his1.class(Error)"). "Error" is a class.
///           And object is name of root object. E.g. his->GetName(). ("his1" in example above).
///         * declaration consist from  properties and values. We try to use the same properties like in CSS where possible
///           in opposite case we use properties from root. E.g. {height: 1200;} for canvases use real css property or
///           {marker-style: 20,21,24,27;} is example of using of property from root;
///   One more important thing - we provide to users opportunity to specify local property
///   in the name of objects, such property has more priority than property from file.
///   E.g. You can define name of object like
///     \code
///         hisArray[i] = new TH1F(TString::Format("his%d.class(Raw).style(marker-color:%d;)", i, i*2 + 3).Data(),
///                           TString::Format("his%d.class(Raw)", i).Data(), 100, -5, 5);
///      \endcode
///    How it made in MakeTestPlot() from $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+.
///
///   The sequance of work with AliDrawStyle:
///     1. Define the style (create file with description of style);
///     2. Register this style;
///     3. Apply it;
///
///   Method of applying the style - AliDrawStyle::ApplyCssStyle() has three input parameters:
///     1. Name of style (it has a type - const char *, so for using TString use .Data());
///     2. TPad object in order to apply your style recursively according with you style;
///     3. Verbosity mode (4 is for detailed output);
///
///   Also we support units and color formats. So, you can use rgb or hex notation for colors in your .css file and pixels or percents for sizes.
///   (See examples)[https://github.com/alisw/AliRoot/blob/master/STAT/test/alirootTestStyle.css]
///  * Usage (work in progress)
///    * performance reports -  with styles as a parameter
///    * QA reports
///    * AliTreePlayer, and TStatToolkit
/// \author [Marian Ivanov](marian.ivanov@cern.ch), [Boris Rumyantsev](boris.rumyantsev@cern.ch)
///
///  ## Example usage
/*!
\code
 .L $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+
 TCanvas *canv = MakeTestPlot(3);
 AliDrawStyle::RegisterCssStyle("test1",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/test1.css",0));
 AliDrawStyle::ApplyCssStyle(canv, "test1");
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
#include "TPad.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TClass.h"
#include "TLegend.h"
#include "TString.h"
#include "TList.h"
#include "TObject.h"
#include "TPad.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TList.h"

#include <iostream>
#include <fstream>
#include <ios>

Int_t AliDrawStyle::padNumber = 0;


/// AliDrawStyle::ParseDeclaration parse input declaration and return values
/// \param input        - input string
/// \param propertyName - name of property to find
/// \return             - property propertyName from the input string using CSS like parsing, empty string in case not found
///
/*!
 ####  Example use:
  TString input="{\nmarker-style:25,21,22,23; \nmarker-color:1,2,4,5; \n}";
  AliDrawStyle::ParseDeclaration(input.Data(),"marker-style");  //  return  "25,21,22,23"
  AliDrawStyle::ParseDeclaration(input.Data(),"marker-color");  //  return  "1,2,4,5"
 */
TString  AliDrawStyle::ParseDeclaration(const char *inputDec, const char *propertyName) {
  TString input(inputDec);
  TPRegexp valPat(":.*;");
  TString dVal = "";
  Int_t index0 = input.Index(propertyName);
  if (index0<0) return "";
  dVal = TString(input(valPat,index0));
  Int_t valLength = dVal.Index(';') - dVal.Index(':') -1;
  return dVal(1,valLength);
}

///  Method returns input number with predefined type
/// \tparam T - type of return value could be any number format from c++
/// \param inputStr - input string could be declaration from css or only string with values
/// \param status - it's a flag which allows to understand should we use the returned value or not
/// \param index - the number of returning value, by default - 0, in case if index more than number of values, last value will be return
/// \param propertyName - name of property, using for sizes
/// \param verbose
/// \param sep - separator from input string. by default - ","
/// \param ignoreBrackets - separators inside of this brackets will be ignore. by default "()"
/// \return - number of type T or -1 if something went wrong
/*!
 ####  Example use:
 \code
 1.
  Bool_t status;
  AliDrawStyle::GetNamedTypeAt<int>("11,13,14", status)
  (int) 11
 2.
  Bool_t status;
  AliDrawStyle::GetNamedTypeAt<int>("11,13,14", status, 1)
  (int) 12
 3.
  Bool_t status;
  AliDrawStyle::GetNamedTypeAt<int>("rgb(12,12,12),#dfdfdf,14", status)
  (int) 1179
 4.
  Bool_t status;
  AliDrawStyle::GetNamedTypeAt<int>("rgb(12,12,12),#dfdfdf,14", status, 1)
  (int) 1180
 5. Bool_t status;
    AliDrawStyle::GetNamedTypeAt<int>("marker-color:rgb(12,12,12),#dfdfdf,14;marker-size:1,2,3,4;", status,2,"marker-size")
    (int) 3
 6. Bool_t status;
    AliDrawStyle::GetNamedTypeAt<int>("marker-color:rgb[12,12,12],#dfdfdf,14;marker-size:1,2,3,4;", status,2,"marker-color", 0, ',',"[]")
    (int) 14
 \endcode
 */
 //TODO: remove template because https://github.com/alisw/AliRoot/pull/657#discussion_r185746846
 Float_t AliDrawStyle::GetNamedTypeAt(const char *inputStr, Bool_t &status, int index, const char *propertyName, Int_t verbose, const char sep, const char *ignoreBrackets) {
  TString inputTStr;
  if(TString(propertyName) != TString("")) inputTStr = AliDrawStyle::ParseDeclaration(inputStr,propertyName);
  else inputTStr = TString(inputStr);
  Float_t res = -1.;
  Int_t arg = 0, startIndex = 0;
  if (TString(inputStr) == TString("")) {
    ::Error("AliDrawStyle", "AliDrawStyle::GetNamedTypeAt(\"%s\", %d, \"%s\"). Options string should not be empty.", inputStr, index, propertyName);
    status = kFALSE;
    return res;
  }
  std::vector<TString> values;
  for (Int_t i = 0; i <= index; i++) values.push_back(TString(""));

  for (Int_t i = 0; i <= inputTStr.Length(); i++) {
    if (arg > index) break;
    if (inputTStr(i) == TString(sep) || i == inputTStr.Length()) {
      values[arg] = TString(inputTStr(startIndex, i - startIndex));
      arg++;
      startIndex = i + 1;
    } else if (inputTStr(i) == TString(ignoreBrackets)[0]) {
      i = inputTStr.Index(ignoreBrackets[1], i);
      continue;
    }
  }

  if (index >= arg && index != 0) return AliDrawStyle::GetNamedTypeAt(inputStr, status, arg-1, propertyName, verbose, sep, ignoreBrackets); // keeping the last index for objects in case index doesn't exist.

  if (values[index].IsDec() || values[index].IsFloat()) {
    if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::GetNamedTypeAt(\"%s\",%d,%d,\"%s\") was transformed in %f",inputStr,status,index,propertyName, values[index].Atof());
    status = kTRUE;
    return values[index].Atof();
  }
  else if (values[index].Contains("%") || values[index].Contains("px")) {
    Float_t valF = AliDrawStyle::ConvertUnit(values[index].Data(),propertyName);
    if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::GetNamedTypeAt(\"%s\",%d,%d,\"%s\") was transformed in %f",inputStr,status,index,propertyName, valF);
    if (valF != -1.0) status = kTRUE;
    return valF;
  }
  else if (values[index].Contains("rgb") || values[index].Contains("#")) {
    Int_t valI = AliDrawStyle::ConvertColor(values[index]);
    if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::GetNamedTypeAt(\"%s\",%d,%d,\"%s\") was transformed in %d",inputStr,status,index,propertyName, valI);
    if (valI != -1) status = kTRUE;
    return (Float_t) valI;
  }
  else {
    status = kFALSE;
    if (verbose == 4) ::Error("AliDrawStyle", "Something went wrong in AliDrawStyle::GetNamedTypeAt(\"%s\",%d,%d,\"%s\"). Result string is \"%s\"",inputStr,status,index,propertyName, values[index].Data());
    return -1.0;
  }

   return -1.0;
}

/// \brief converter from pixels to float
///  Pixels are relative units, it means that we should have something for converting from it.
///   * In case margins we use gPad->GetWw() for left-right margin, and gPad->GetWh() for top-bottom margins.
///   * In case markers-size we use this:
///       The marker size does not refer to any coordinate systems, it is an absolute value.
///       Therefore the marker size is not affected by any change in TPad's scale.
///       A marker size equl to 1 correspond to 8 pixels. That is, a square marker with size 1
///       will be drawn with a side equal to 8 pixels on the screen.
///       https://root.cern.ch/root/html600/TAttMarker.html
/// \param value - value for converting
/// \param option - sets way of the converting
/// \param verbose
/// \return - Float_t number or -1.0 if something went wrong
Float_t AliDrawStyle::PixelsToFloat_t(const char *value, const char *option, Int_t verbose) {
  TString valStr = TString(value);
  TString optionT = TString(option);
  Float_t val = (Float_t) -1.0;

  if (optionT.Contains("marker")) {
    val = (Float_t) TString(valStr(0, valStr.Length() - 2)).Atof() / 8;
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::PixelsToFloat_t(\"%s\",\"%s\") transformed into %f.",value, option,val);
    return val;
  }

  if (gPad == nullptr) {
    ::Error("AliDrawStyle", "AliDrawStyle::PixelsToFloat_t(\"%s\",\"%s\"). For converting of pixels you must have a TPad object.", value, option);
    return (Float_t) -1.0;
  }
  else {
    if (optionT.Contains("top") || optionT.Contains("bottom")) {
      val = (Float_t) TString(valStr(0, valStr.Length() - 2)).Atof() / gPad->GetWh();
    } else if (optionT.Contains("left") || optionT.Contains("right")) {
      val = (Float_t) TString(valStr(0, valStr.Length() - 2)).Atof() / gPad->GetWw();
    } else if (optionT.Contains("border")) {
      val = (Float_t) TString(valStr(0, valStr.Length() - 2)).Atof();
    }
    if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::PixelsToFloat_t(\"%s\",\"%s\") transformed into %f. Calculations based on pad - %s", value, option, val, gPad->GetName());
  }

  return val;
}

/// \brief converter from percents to float
/// \param value - value for converting
/// \param verbose
/// \return - Float_t number or -1.0 if something went wrong
Float_t AliDrawStyle::PercentToFloat_t(const char *value, Int_t verbose) {
  TString valStr = TString(value);
  Float_t val = (Float_t) TString(valStr(0, valStr.Length() - 1)).Atof() / 100;
  if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::PercentToFloat_t(%s) transformed into %f", valStr.Data(), val);
  return val;
}

/// \brief Defines what units user used and call appropriate converter
/// Rules for values:
///   *for pixels we are using "px" in the end of value, so if you want to recieve 300 pixels, just set "300px"
///    as input value
///   * for percents use %. "30%"
/// \param value - value for convert
/// \param option
/// \return \return - Float_t number or -1.0 if something went wrong
Float_t AliDrawStyle::ConvertUnit(const char *inputValues, const char *option, Int_t verbose) {
  TString value = TString(inputValues);
  if (value.Contains("px") && option != "")
    return AliDrawStyle::PixelsToFloat_t(value, option, verbose);
  else if (value.Contains("%"))
    return AliDrawStyle::PercentToFloat_t(value, verbose);
  else if (value.IsFloat()) return (Float_t) value.Atof();
  if (verbose == 4) ::Error("AliDrawStyle", "In AliDrawStyle::ConvertUnit(%s,%s) occured the error.", inputValues, option);
  return (Float_t) -1.0;
}

/// \brief converter from RGB to Int_t (root format of colors)
/// \param inputColor
/// \param verbose
/// \return
/*!
 ####  Example use:
 \code
    AliDrawStyle::RgbToColor_t("rgb(0,0,0)")
    (int) 1
    AliDrawStyle::RgbToColor_t("rgb(255,255,255)")
    (int) 0
 \endcode
 */
Int_t AliDrawStyle::RgbToColor_t(const char *inputString, Int_t verbose) {
  TString rgbValuesStr = "";
  TPRegexp rgbPattern("[(].*[)]");
  TString inputColor = TString(inputString);
  if (inputColor.CountChar(',') != 2) {
    ::Error("AliDrawStyle", "AliDrawStyle::RgbToColor_t(%s) - rgb should be rgb(Float_t, Float_t, Float_t)",
            inputColor.Data());
    return -1;
  }
  rgbValuesStr = TString(inputColor(rgbPattern)).ReplaceAll("(", "");
  rgbValuesStr = rgbValuesStr.ReplaceAll(")", "");
  TObjArray *rgbValues = rgbValuesStr.Tokenize(",");
  if (TString(rgbValues->At(0)->GetName()).IsDec() && TString(rgbValues->At(1)->GetName()).IsDec() &&
      TString(rgbValues->At(2)->GetName()).IsDec()) {
    Int_t r = TString(rgbValues->At(0)->GetName()).Atoi();
    Int_t g = TString(rgbValues->At(1)->GetName()).Atoi();
    Int_t b = TString(rgbValues->At(2)->GetName()).Atoi();
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::RgbToColor_t(%s) transformed into "
                     "TColor::GetColor(%d, %d, %d) and returned %d", inputColor.Data(), r, g, b,
             TColor::GetColor(r, g, b));
    return TColor::GetColor(r, g, b);
  }
  else if (TString(rgbValues->At(0)->GetName()).IsFloat() && TString(rgbValues->At(1)->GetName()).IsFloat() &&
           TString(rgbValues->At(2)->GetName()).IsFloat()) {
    Float_t r = (Float_t) TString(rgbValues->At(0)->GetName()).Atof();
    Float_t g = (Float_t) TString(rgbValues->At(1)->GetName()).Atof();
    Float_t b = (Float_t) TString(rgbValues->At(2)->GetName()).Atof();
    if (verbose == 4) ::Info("AliDrawStyle","AliDrawStyle::RgbToColor_t(%s) transformed into "
                                     "TColor::GetColor(%f, %f, %f) and returned %d", inputColor.Data(),r,g,b,
                             TColor::GetColor(r, g, b));
    return TColor::GetColor(r, g, b);
  }
  if (verbose == 4)
    ::Warning("AliDrawStyle", "In AliDrawStyle::RgbToColor_t(%s) something went wrong", inputString);
  return -1;
}

/// \brief converter from HEX to Int_t (root format of colors)
/// \param inputColor
/// \param verbose
/// \return
Int_t AliDrawStyle::HexToColor_t(const char *inputString, Int_t verbose) {
  TString inputColor = TString(inputString);
  if (verbose == 4) ::Info("AliDrawStyle","AliDrawStyle::HexToColor_t(%s) transformed into "
                                   "TColor::GetColor(%s) and returned %d", inputColor.Data(), inputColor.Data(),
                           TColor::GetColor(inputColor.Data()));
  return TColor::GetColor(inputColor.Data());
}

/// \brief Defines what format of color user used and call appropriate converter
/// Rules for values:
///   *for rgb we are using "rgb()" in the end of value, so if you want to recieve red color, just set "rgb(255,0,0)"
///    as input value
///   * for hex use "#ff0000"
/// \param inputString - string of input values ("#000000", "rgb(0,0,0)", "0")
/// \param verbose
/// \return - TColor::GetColor(value)
/*!
 ####  Example use:
 \code
    AliDrawStyle::HexToColor_t("#ff0000")
    (int) 2
    AliDrawStyle::HexToColor_t("#00ff00")
    (int) 3
 \endcode
 */
Int_t AliDrawStyle::ConvertColor(const char *inputString, Int_t verbose) {
  TString value = TString(inputString);
  Int_t color;
  if (value(0,3) == TString("rgb"))
    color = AliDrawStyle::RgbToColor_t(value, verbose);
  else if (value(0) == TString("#"))
    color = AliDrawStyle::HexToColor_t(value, verbose);
  else if (value.IsDec())
    color = value.Atoi();
  else {
    ::Error("AliDrawStyle", "AliDrawStyle::ConvertColor(%s) - something wrong with colors", value.Data());
    return -1;
  }
  if (verbose == 4) ::Info("AliDrawStyle","AliDrawStyle::ConvertColor(%s) transformed into "
            "returned %d", inputString,color);
  return color;
}

/// \brief prepares value for applying.
///  This method defines which value user specified (color, unit, number, etc) and return number with appropriate type
///  Returning value will use in Set methods of objects. (TH1.SetMarkerColor(value)).
/// \tparam T
/// \param styleName - name of preloading style
/// \param propertyName - name of property which you want to change
/// \param elementID - root class of object
/// \param classID  - tag or user predefined class
/// \param objectID - name of root-object
/// \param localStyle - local style from name of object with high priority
/// \param objNum - number of object in pad
/// \param verbose
/// \return
//TODO: does it look like godness function? refactor? @Boris
//TODO: perhaps I should combine {TString elementID, TString classID, TString objectID, TString localStyle} to the array? @Boris
Float_t AliDrawStyle::PrepareValue(const char* styleName, TString propertyName, TString elementID, TString classID, TString objectID, TString localStyle, Bool_t &status, Int_t objNum, Int_t verbose) {
  TString property = "";
  TString cProperty = "";
  Int_t hisNum = -1;
  status = kFALSE;
  Float_t value;
  if (localStyle.Contains(propertyName.Data()) && (propertyName.Contains("size") || propertyName.Contains("margin"))) {
    property = localStyle;
    hisNum = 0;
    cProperty = propertyName;
    if (verbose == 4) ::Info("AliDrawStyle", "1.AliDrawStyle::PrepareValue(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d) provided value = \"%s\" ", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data(), localStyle.Data(), objNum, property.Data());
  }
  else if (localStyle.Contains(propertyName.Data())) {
    property = AliDrawStyle::ParseDeclaration(localStyle.Data(), propertyName.Data());
    hisNum = 0;
    cProperty = "";
    if (verbose == 4) ::Info("AliDrawStyle", "2.AliDrawStyle::PrepareValue(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d) provided value = \"%s\" ", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data(), localStyle.Data(), objNum, property.Data());
  }
  else if (propertyName.Contains("size") || propertyName.Contains("margin")) {
    property = AliDrawStyle::GetValue(styleName, "", elementID, classID, objectID, "", verbose);
    hisNum = objNum;
    cProperty = propertyName;
    if (!property.Contains(cProperty)) {
      if (verbose == 4) ::Info("AliDrawStyle","AliDrawStyle::PrepareValue(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\") property not found in css file", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data());
      return -1.;
    }
    if (verbose == 4) ::Info("AliDrawStyle", "3.AliDrawStyle::PrepareValue(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d) provided value = \"%s\" ", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data(), localStyle.Data(), objNum, property.Data());
  }
  else {
    property = AliDrawStyle::GetValue(styleName, propertyName, elementID, classID, objectID, localStyle, verbose);
    hisNum = objNum;
    cProperty = "";
    if (property == TString("")) return -1.;
    if (verbose == 4) ::Info("AliDrawStyle", "4.AliDrawStyle::PrepareValue(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d) provided value = \"%s\" ", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data(), localStyle.Data(), objNum, property.Data());
  }

  value = AliDrawStyle::GetNamedTypeAt(property.Data(), status, hisNum, cProperty, verbose);
  if (verbose == 4) ::Info("AliDrawStyle", "5.AliDrawStyle::PrepareValue(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d) status is %d and returned value", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data(), localStyle.Data(), objNum, status);
  if (status) return value;

  return -1.0;
}

/// Read CSS html like files  (*see also AliRoot modification in CSS)
/// TODO:
/// * proper exception  handling (Boris)
///   * Use ::Error verbosity according debug level
/// \param inputName     - input file to read
/// \param verbose       - specify verbose level for ::error and ::info (Int_t should be interpreted as an bit-mask)
/// \return              - TObjArray  with the pairs TNamed of the CSS <Selector, declaration> or  TObjArray (recursive structure like includes)
TObjArray *AliDrawStyle::ReadCSSFile(const char *  inputName, TObjArray * cssArray, Int_t verbose) {
  //check file exist
  if (gSystem->GetFromPipe(TString("[ -f ") + TString(inputName) +  TString(" ] && echo 1 || echo 0")) == "0") {
    ::Error("AliDrawStyle::ReadCSSFile","File %s doesn't exist", inputName);
    return nullptr;
  }
  TString inputCSS = gSystem->GetFromPipe(TString::Format("cat %s",inputName).Data());     // I expect this variable is defined
  //remove comments:
  while (inputCSS.Index("*/") > 0) {
    inputCSS = inputCSS(0, inputCSS.Index("/*")) + inputCSS(inputCSS.Index("*/")+2, inputCSS.Length());
  }
  //inputCSS.ReplaceAll("\n", ""); //  check performance and implement better variant;
  TObjArray *tokenArray = inputCSS.Tokenize("{}");   //assuming we can not use {} symbols in the style IDS
  Int_t entries = tokenArray->GetEntries();
  if (cssArray==nullptr) {
    cssArray = new TObjArray(entries / 2);
  }
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
/// \param pCssOut     - output stream ( )
void AliDrawStyle::WriteCSSFile(TObjArray * cssArray, const char *  outputName, std::fstream *pCssOut) {

  if (pCssOut == nullptr) {
    pCssOut=new std::fstream;
    pCssOut->open(outputName, std::ios_base::out|std::ios_base::trunc);
  }
  std::fstream &cssOut = *pCssOut;
  for (Int_t i=0;i<cssArray->GetEntries();i++) {
    TObject *object = cssArray->At(i);
    if (object->InheritsFrom("TObjArray")) {
      AliDrawStyle::WriteCSSFile((TObjArray*)object, nullptr, pCssOut);
    }else {
      cssOut << object->GetName();
      cssOut << "{";
      cssOut << object->GetTitle();
      cssOut << "}";
    }
  }
  cssOut<<std::flush;
  if (outputName!=nullptr) {
    pCssOut->close();
    delete pCssOut;
  }
}
//TODO: combine ElementSearch, ClassSearch, ObjectSearch into one function
/// \brief Checks element id in selector
/// \param selector - input string with selectors
/// \param elementID
/// \param verbose
/// \return true if element was found
Bool_t AliDrawStyle::ElementSearch(const TString selector, const TString elementID, Int_t verbose) {
  if (selector == TString("")) {
    ::Error("AliDrawStyle", "AliDrawStyle::ElementSearch(\"%s\", \"%s\") selector should not be empty", selector.Data(), elementID.Data());
    return kFALSE;
  }
  //impossible for us
  if (elementID == TString("")) {
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::ElementSearch(\"%s\", \"%s\") elementID is empty, any element is fine", selector.Data(), elementID.Data());
    return kTRUE;
  }
  TString elementFromSelector = "";
  Int_t finish = -1;
  if(selector.Index('#') >= 0) finish = selector.Index('#');
  if(selector.Index('.') >= 0) finish = selector.Index('.');

  if (finish == 0) return kTRUE;
  else if (finish > 0) elementFromSelector = selector(0, finish);
  else elementFromSelector = selector;
  if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::ElementSearch(\"%s\", \"%s\"). Selector was transformed to %s", selector.Data(), elementID.Data(), elementFromSelector.Data());

  if (elementFromSelector.Index('*') >= 0) {
    elementFromSelector = elementFromSelector.ReplaceAll("*", ".*");
    TPRegexp elemPattern(elementFromSelector);
    if (elementID(elemPattern) == elementID || TClass(elementID.Data()).InheritsFrom(TString(elementFromSelector(0,elementFromSelector.Index("."))).Data())) {
      if (verbose == 4)
        ::Info("AliDrawStyle", "AliDrawStyle::ElementSearch(\"%s\", \"%s\") was successful", selector.Data(), elementID.Data());
      return kTRUE;
    }
    else {
      if (verbose == 4)
        ::Info("AliDrawStyle", "AliDrawStyle::ElementSearch(\"%s\", \"%s\") was unsuccessful", selector.Data(), elementID.Data());
      return kFALSE;
    }
  }

  if (elementFromSelector == elementID) {
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::ElementSearch(\"%s\", \"%s\") was successful", selector.Data(), elementID.Data());
    return kTRUE;
  }
  if (verbose == 4)
    ::Info("AliDrawStyle", "AliDrawStyle::ElementSearch(\"%s\", \"%s\") was unsuccessful", selector.Data(), elementID.Data());
  return kFALSE;
}

/// \brief Checks class in selector
/// \param selector - input string with selectors
/// \param classID
/// \param verbose
/// \return true if class was found
Bool_t AliDrawStyle::ClassSearch(const TString selector, const TString classID, Int_t verbose) {
  if (selector == TString("")) {
    ::Error("AliDrawStyle", "AliDrawStyle::ClassSearch(\"%s\", \"%s\") selector should not be empty", selector.Data(), classID.Data());
    return kFALSE;
  }
  if (classID == TString("")) {
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::ClassSearch(\"%s\", \"%s\") classID is empty, any class is fine", selector.Data(), classID.Data());
    return kTRUE;
  }

  TPRegexp classPattern("[.].*#|[.].*");
  TString classFromSelector = "";

  classFromSelector = selector(classPattern);
  if (classFromSelector == TString(""))
    return kTRUE;

  if(classFromSelector.Index('#') >= 0) classFromSelector = classFromSelector(1,classFromSelector.Length() - 2);
  else classFromSelector = classFromSelector(1,classFromSelector.Length());
  if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::ClassSearch(\"%s\", \"%s\"). Selector was transformed to %s", selector.Data(), classID.Data(), classFromSelector.Data());

  //for multi classes
  TObjArray *classIDs = classID.Tokenize(",");
  Int_t nC = classIDs->GetEntriesFast();;
  TString tempClassID = "";

  for (Int_t i = 0; i < nC; i++) {
    tempClassID = classIDs->At(i)->GetName();
    if (classFromSelector.Index('*') >= 0) {
      classFromSelector = classFromSelector.ReplaceAll("*", ".*");
      TPRegexp classSelPattern(classFromSelector);
      if (tempClassID(classSelPattern) == tempClassID) {
        if (verbose == 4)
          ::Info("AliDrawStyle", "AliDrawStyle::ClassSearch(\"%s\", \"%s\") was successful", selector.Data(), classID.Data());
        return kTRUE;
      }
      else {
        ::Info("AliDrawStyle", "AliDrawStyle::ClassSearch(\"%s\", \"%s\") was unsuccessful", selector.Data(), classID.Data());
        return kFALSE;
      }
    }

    if (classFromSelector == tempClassID) {
      if (verbose == 4)
        ::Info("AliDrawStyle", "AliDrawStyle::ClassSearch(\"%s\", \"%s\") was successful", selector.Data(), classID.Data());
      return kTRUE;
    }
  }
  if (verbose == 4)
    ::Info("AliDrawStyle", "AliDrawStyle::ClassSearch(\"%s\", \"%s\") was unsuccessful", selector.Data(), classID.Data());
  return kFALSE;
}

///  Checks object in selector
/// \param selector
/// \param objectID
/// \param verbose
/// \return true if object was found
Bool_t AliDrawStyle::ObjectSearch(const TString selector, const TString objectID, Int_t verbose) {
  if (selector == TString("")) {
    ::Error("AliDrawStyle", "AliDrawStyle::ObjectSearch(\"%s\", \"%s\") selector should not be empty", selector.Data(), objectID.Data());
    return kFALSE;
  }
  if (objectID == TString("")) {
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::ObjectSearch(\"%s\", \"%s\") objectID is empty.", selector.Data(), objectID.Data());
    return kTRUE;
  }

  TPRegexp objectPattern("[#].*");
  TString objectFromSelector = "";

  objectFromSelector = selector(objectPattern);
  if (objectFromSelector == TString(""))
    return kTRUE;

  objectFromSelector = objectFromSelector(1,objectFromSelector.Length());
  if (verbose == 4) ::Info("AliDrawStyle", "AliDrawStyle::ObjectSearch(\"%s\", \"%s\"). Selector was transformed to %s", selector.Data(), objectID.Data(), objectFromSelector.Data());

  if (objectFromSelector.Index('*') >= 0) {
    objectFromSelector = objectFromSelector.ReplaceAll("*", ".*");
    TPRegexp objectSelPattern(objectFromSelector);
    if (objectID(objectSelPattern) == objectID) {
      if (verbose == 4)
        ::Info("AliDrawStyle", "AliDrawStyle::ObjectSearch(\"%s\", \"%s\") was successful", selector.Data(), objectID.Data());
      return kTRUE;
    }
    else {
      if (verbose == 4)
        ::Info("AliDrawStyle", "AliDrawStyle::ObjectSearch(\"%s\", \"%s\") was unsuccessful", selector.Data(), objectID.Data());
      return kFALSE;
    }
  }

  if (objectFromSelector == objectID) {
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::ObjectSearch(\"%s\", \"%s\") was successful", selector.Data(), objectID.Data());
    return kTRUE;
  }
  if (verbose == 4)
    ::Info("AliDrawStyle", "AliDrawStyle::ObjectSearch(\"%s\", \"%s\") was unsuccessful", selector.Data(), objectID.Data());
  return kFALSE;
}

/// Function of checking the match between "CSS" selector and object parameters: elementID, classID, objectID
/// \param selectors   - selector ID
/// \param elementName - name of element
/// \param className   - name of class
/// \param objectName  - object name
/// \return            - kTRUE if selector match class name and object name
/*!
 ####  Example use:
 \code
  TString selectors = "\n\n\nTH1.Status#obj1, TH1.Warning#obj1, TH1.Warning#obj3 \tTGraph#obj1, TGraph.Status#TPC.QA.dcar_posA_1 \tTGraph.Warning#TPC.QA.dcar_posA_2 \tTF1.Status, .Status#obj1, #obj3";
  AliDrawStyle::IsSelected(selectors, "", "", "obj3")            // return true
  AliDrawStyle::IsSelected(selectors, "TF1", "Warning", "obj3")  // return true
  AliDrawStyle::IsSelected(selectors, "TH1C", "Warning", "obj1") // return false
  \endcode
 */
Bool_t  AliDrawStyle::IsSelected(TString selectors, const TString elementID, const TString classID, const TString objectID, Int_t verbose) {
  if (elementID == TString("") && classID == TString("") && objectID == TString("")) {
    if (verbose == 4)
      ::Info("AliDrawStyle", "AliDrawStyle::IsSelected(\"%s\", \"%s\", \"%s\", \"%s\") returned kFALSE. All elements are empty.", \
                  selectors.Data(), elementID.Data(), classID.Data(), objectID.Data());
    return kFALSE;
  }
  Ssiz_t    fromStart1 = 0;
  Ssiz_t    fromStart     = 0;
  TString   subSelectors   = "";
  TString   selector       = "";
  selectors = selectors.ReplaceAll("\n", "");
  while (selectors.Tokenize(subSelectors, fromStart1, " \t")) {
    fromStart = 0;
    while (subSelectors.Tokenize(selector, fromStart, ", ")) {
      selector = selector.Strip(TString::kBoth, ' ');
      if (ElementSearch(selector, elementID, verbose) && ClassSearch(selector, classID, verbose) && ObjectSearch(selector, objectID, verbose)) {
        if (verbose == 4)
          ::Info("AliDrawStyle", "AliDrawStyle::IsSelected(\"%s\", \"%s\", \"%s\", \"%s\") returned kTRUE", \
                  selectors.Data(), elementID.Data(), classID.Data(), objectID.Data());
        return kTRUE;
      }
    }
  }
  if (verbose == 4)
    ::Info("AliDrawStyle", "AliDrawStyle::IsSelected(\"%s\", \"%s\", \"%s\", \"%s\") returned kFALSE", \
            selectors.Data(), elementID.Data(), classID.Data(), objectID.Data());
  return kFALSE;
}

/// \brief GetValue gets value from css file or local style.
///
/// \param styleName      - name of predefined array with styles (see. AliDrawStyle::RegisterCssStyle())
/// \param propertyName   - name of property according to css notation
/// \param elementID
/// \param classID
/// \param objectID
/// \return TString with values according to input propertyName
/*!
  ####  Example use:
   \code
   AliDrawStyle::RegisterCssStyle("alirootTestStyle.css",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css",0));
   AliDrawStyle::GetValue("alirootTestStyle.css","marker-size", "TGraph", "", "obj1");                      // "33,34,35,36"
   AliDrawStyle::GetValue("alirootTestStyle.css","marker-size", "TGraph", "Status", "TPC.QA.dcar_posA_1");        // "1,2,3,4"
   AliDrawStyle::GetValue("alirootTestStyle.css","marker-size", "TGraph", "Warning", "TPC.QA.dcar_posA_1"); // ""
    \endcode
*/
TString AliDrawStyle::GetValue(const char *styleName, TString propertyName, TString elementID, TString classID, TString objectID, TString localStyle, Int_t verbose) {
  TString value       = "";
  if (localStyle.Contains(propertyName.Data()) && propertyName != TString("")) {
    value = AliDrawStyle::ParseDeclaration(localStyle.Data(), propertyName.Data());
    if (verbose == 4) ::Info("AliDrawStyle","AliDrawStyle::GetValue(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\") was found in local property and was transformed into \"%s\"", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data(), value.Data());
    return value;
  }

  if (fCssStyleAlice[styleName] == nullptr || (elementID == TString("") && classID == TString("") && objectID == TString(""))) return "";

  TString actProperty = "";
  Int_t   entries     = fCssStyleAlice[styleName]->GetEntriesFast();
  TString declaration = "";

  for(Int_t i = 0; i < entries; i++) {
    if (AliDrawStyle::IsSelected(TString(fCssStyleAlice[styleName]->At(i)->GetName()),
                                 elementID, classID, objectID, verbose)) {
      if (propertyName == "") value = fCssStyleAlice[styleName]->At(i)->GetTitle();
      else value = AliDrawStyle::ParseDeclaration(fCssStyleAlice[styleName]->At(i)->GetTitle(), propertyName.Data());
      if (value != "") actProperty = value.Strip(TString::kBoth, ' ');
    }
  }
  if (verbose == 4 && propertyName != TString("")) {
    if (value != TString("")) ::Info("AliDrawStyle","AliDrawStyle::GetValue(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\") was transformed into \"%s\"", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data(), actProperty.Data());
    else ::Info("AliDrawStyle","AliDrawStyle::GetValue(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\") property not found in css file", styleName, propertyName.Data(), elementID.Data(), classID.Data(), objectID.Data());
  }
  return actProperty;
}

///
/// \param styleName    - name of predefined array with styles (see. AliDrawStyle::RegisterCssStyle())
/// \param cGraph       - object inherited from TGraph in which you want to change style
/*!
  ####  Example use:
  See AliDrawStyle::ApplyCssStyle();
*/
void AliDrawStyle::TGraphApplyStyle(const char* styleName, TGraph *cGraph, Int_t objNum, Int_t verbose) {

  AliDrawStyle::ObjectApplyStyle(styleName, cGraph, objNum, verbose);

  Bool_t status = kFALSE;
  TString elementID = "";
  TString classID = "";
  TString objectID = "";
  TString property = "";
  TString localStyle = "";
  Int_t valueI = 0;
  Float_t valueF;

  AliDrawStyle::GetIds(cGraph, elementID, classID, objectID, localStyle, verbose);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("axis-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetAxisColor((Color_t) valueI);
    cGraph->GetYaxis()->SetAxisColor((Color_t) valueI);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("label-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetLabelColor((Color_t) valueI);
    cGraph->GetYaxis()->SetLabelColor((Color_t) valueI);
  }
  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("label-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetLabelFont((Style_t) valueI);
    cGraph->GetYaxis()->SetLabelFont((Style_t) valueI);
  }

  valueF = AliDrawStyle::PrepareValue(styleName, TString("label-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetLabelSize(valueF);
    cGraph->GetYaxis()->SetLabelSize(valueF);
  }

  valueF = AliDrawStyle::PrepareValue(styleName, TString("label-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetLabelOffset(valueF);
    cGraph->GetYaxis()->SetLabelOffset(valueF);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("title-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetTitleFont((Style_t) valueI);
    cGraph->GetYaxis()->SetTitleFont((Style_t) valueI);
  }

  valueF = AliDrawStyle::PrepareValue(styleName, TString("title-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetTitleOffset(valueF);
    cGraph->GetYaxis()->SetTitleOffset(valueF);
  }

  valueF = AliDrawStyle::PrepareValue(styleName, TString("title-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cGraph->GetXaxis()->SetTitleSize(valueF);
    cGraph->GetYaxis()->SetTitleSize(valueF);
  }
}

///
/// \param styleName    - name of predefined array with styles (see. AliDrawStyle::RegisterCssStyle())
/// \param cHis       - object inherited from TH1 in which you want to change style
/*!
  ####  Example use:
  See AliDrawStyle::ApplyCssStyle();
*/
void AliDrawStyle::TH1ApplyStyle(const char *styleName, TH1 *cHis, Int_t objNum, Int_t verbose) {

  AliDrawStyle::ObjectApplyStyle(styleName, cHis, objNum, verbose);

  Bool_t status = kFALSE;
  TString elementID = "";
  TString classID = "";
  TString objectID = "";
  TString property = "";
  TString localStyle = "";
  Int_t valueI = 0;
  Float_t valueF;

  AliDrawStyle::GetIds(cHis, elementID, classID, objectID, localStyle, verbose);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("axis-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetAxisColor((Color_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("label-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelColor((Color_t) valueI, "xyz");

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("Xlabel-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelColor((Color_t) valueI, "X");

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("Ylabel-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelColor((Color_t) valueI, "Y");

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("Zlabel-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelColor((Color_t) valueI, "Y");

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("label-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelFont((Style_t) valueI, "xyz");

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("Xlabel-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelFont((Style_t) valueI, "X");

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("Ylabel-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelFont((Style_t) valueI, "Y");

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("Zlabel-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelFont((Style_t) valueI, "Z");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("label-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelSize(valueF, "xyz");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Xlabel-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelSize(valueF, "X");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Ylabel-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelSize(valueF, "Y");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Zlabel-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelSize(valueF, "Z");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("label-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetLabelOffset(valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("title-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleFont((Style_t) valueI);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("title-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleOffset(valueF, "xyz");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Xtitle-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleOffset(valueF, "X");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Ytitle-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleOffset(valueF, "Y");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Ztitle-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleOffset(valueF, "Z");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("title-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleSize(valueF, "xyz");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Xtitle-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleSize(valueF, "X");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Ytitle-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleSize(valueF, "Y");

  valueF = AliDrawStyle::PrepareValue(styleName, TString("Ztitle-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cHis->SetTitleSize(valueF, "Z");

  //TODO: how to make it not only for TF1? @Boris
  for (Int_t i = 0; i < cHis->GetListOfFunctions()->GetEntries(); i++) {
    AliDrawStyle::ObjectApplyStyle(styleName, (TF1 *) cHis->GetListOfFunctions()->At(i), objNum, verbose);
  }

}

///
/// \param styleName    - name of predefined array with styles (see. AliDrawStyle::RegisterCssStyle())
/// \param cFunc       - object inherited from TF1 in which you want to change style
/*!
  ####  Example use:
  See AliDrawStyle::ApplyCssStyle();
  //still not implemented applyStyle for fitFunction for TH1
*/
void AliDrawStyle::TF1ApplyStyle(const char *styleName, TF1 *cFunc, Int_t objNum, Int_t verbose) {

  AliDrawStyle::ObjectApplyStyle(styleName, cFunc, objNum, verbose);

  Bool_t status = false;
  TString elementID = "";
  TString classID = "";
  TString objectID = "";
  TString property = "";
  TString localStyle = "";
  Int_t valueI = 0;

  AliDrawStyle::GetIds(cFunc, elementID, classID, objectID, localStyle, verbose);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("axis-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cFunc->GetXaxis()->SetAxisColor((Color_t) valueI);
    cFunc->GetYaxis()->SetAxisColor((Color_t) valueI);
    if (cFunc->GetZaxis()) cFunc->GetZaxis()->SetAxisColor((Color_t) valueI);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("label-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cFunc->GetXaxis()->SetLabelFont((Style_t) valueI);
    cFunc->GetYaxis()->SetLabelFont((Style_t) valueI);
    if (cFunc->GetZaxis()) cFunc->GetZaxis()->SetLabelFont((Style_t) valueI);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("label-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cFunc->GetXaxis()->SetLabelSize(valueI);
    cFunc->GetYaxis()->SetLabelSize(valueI);
    if (cFunc->GetZaxis()) cFunc->GetZaxis()->SetLabelSize(valueI);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("label-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cFunc->GetXaxis()->SetLabelOffset(valueI);
    cFunc->GetYaxis()->SetLabelOffset(valueI);
    if (cFunc->GetZaxis()) cFunc->GetZaxis()->SetLabelOffset(valueI);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("title-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cFunc->GetXaxis()->SetTitleFont((Style_t) valueI);
    cFunc->GetYaxis()->SetTitleFont((Style_t) valueI);
    if (cFunc->GetZaxis()) cFunc->GetZaxis()->SetTitleFont((Style_t) valueI);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("title-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cFunc->GetXaxis()->SetTitleSize(valueI);
    cFunc->GetYaxis()->SetTitleSize(valueI);
    if (cFunc->GetZaxis()) cFunc->GetZaxis()->SetTitleSize(valueI);
  }

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("title-offset"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) {
    cFunc->GetXaxis()->SetTitleOffset(valueI);
    cFunc->GetYaxis()->SetTitleOffset(valueI);
    if (cFunc->GetZaxis()) cFunc->GetZaxis()->SetTitleOffset(valueI);
  }
}

///
/// \param styleName
/// \param cLegend
/// \param objNum
/// \param verbose
void AliDrawStyle::TLegendApplyStyle(const char *styleName, TLegend *cLegend, Int_t objNum, Int_t verbose) {
  Bool_t status = false;
  TString elementID = "";
  TString classID = "";
  TString objectID = "";
  TString property = "";
  TString localStyle = "";
  Int_t valueI;
  Float_t valueF;
  Double_t valueD;

  AliDrawStyle::GetIds(cLegend, elementID, classID, objectID, localStyle, verbose);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("column-separation"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetColumnSeparation(valueF);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetMargin(valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("ncolumns"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetNColumns(valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("border-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetBorderSize(valueI);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("corner-radius"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetCornerRadius(valueD);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("shadow-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetShadowColor(valueI);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("x1-ndc"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetX1NDC(valueD);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("x2-ndc"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetX2NDC(valueD);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("y1-ndc"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetY1NDC(valueD);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("y2-ndc"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetY2NDC(valueD);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("x1"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetX1(valueD);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("x2"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetX2(valueD);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("y1"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetY1(valueD);

  valueD = (Double_t) AliDrawStyle::PrepareValue(styleName, TString("y2"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetY2(valueD);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("bbox-centerx"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetBBoxCenterX(valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("bbox-centery"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetBBoxCenterY(valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("bbox-x1"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetBBoxX1(valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("bbox-x2"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetBBoxX2(valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("bbox-y1"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetBBoxY1(valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("bbox-y2"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetBBoxY2(valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("line-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetLineColor((Color_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("line-style"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetLineStyle((Style_t) valueI);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("line-width"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetLineWidth((Width_t) valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("fill-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetFillColor((Color_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("fill-style"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetFillStyle((Style_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("text-align"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetTextAlign((Short_t) valueI);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("text-angle"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetTextAngle(valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("text-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetTextColor((Color_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("text-font"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetTextFont((Font_t) valueI);


  valueF = AliDrawStyle::PrepareValue(styleName, TString("text-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cLegend->SetTextSize(valueF);


/*
SetEntryLabel
SetEntryOption
SetEntrySeparation
SetHeader
SetToolTipText
SetDrawOption
SetDtorOnly
SetObjectStat

*SetTextSizePixels
 */
}

/// \brief Common part of TGraphApplyStyle(), TH1ApplyStyle(), TF1ApplyStyle()
/// \tparam T
/// \param styleName
/// \param cObj
/// \param objNum
/// \param verbose
/// \return
template <typename T>
 void AliDrawStyle::ObjectApplyStyle(const char *styleName, T *cObj, Int_t objNum, Int_t verbose) {

  Bool_t status = false;
  TString elementID = "";
  TString classID = "";
  TString objectID = "";
  TString property = "";
  TString localStyle = "";
  Int_t valueI;
  Float_t valueF;

  AliDrawStyle::GetIds(cObj, elementID, classID, objectID, localStyle, verbose);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("marker-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetMarkerColor((Color_t) valueI);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("marker-size"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetMarkerSize(valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("marker-style"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetMarkerStyle((Style_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("line-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetLineColor((Color_t) valueI);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("line-width"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetLineWidth((Width_t) valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("line-style"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetLineStyle((Style_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("fill-color"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetFillColor((Color_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("fill-style"), elementID, classID, objectID, localStyle, status, objNum, verbose);
  if (status) cObj->SetFillStyle((Style_t) valueI);
}

///
/// \param styleName    - name of predefined array with styles (see. AliDrawStyle::RegisterCssStyle())
/// \param cPad       - object inherited from TPad in which you want to change style
/*!
  ####  Example use:
  See AliDrawStyle::ApplyCssStyle();
*/
void AliDrawStyle::TPadApplyStyle(const char *styleName, TPad *cPad, Int_t verbose) {

  Bool_t status = false;
  TString elementID = "";
  TString classID = "";
  TString objectID = "";
  TString property = "";
  TString localStyle = "";
  Int_t valueI;
  Float_t valueF;

  AliDrawStyle::GetIds(cPad, elementID, classID, objectID, localStyle, verbose);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("fill-color"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetFillColor((Color_t) valueI);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) {
    cPad->SetBottomMargin(valueF);
    cPad->SetTopMargin(valueF);
    cPad->SetLeftMargin(valueF);
    cPad->SetRightMargin(valueF);
  }
  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-bottom"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetBottomMargin(valueF);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-top"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetTopMargin(valueF);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-left"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetLeftMargin(valueF);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-right"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetRightMargin(valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("border-size"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetBorderSize((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("border-mode"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetBorderMode((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("gridX"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetGridx((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("gridY"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetGridy((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("tickX"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetTickx((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("tickY"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetTicky((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("logX"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetLogx((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("logY"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetLogy((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("logZ"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cPad->SetLogz((Short_t) valueI);
}

///
/// \param styleName    - name of predefined array with styles (see. AliDrawStyle::RegisterCssStyle())
/// \param cCanvas       - object inherited from TCanvas in which you want to change style
/*!
  ####  Example use:
  See AliDrawStyle::ApplyCssStyle();
*/
void AliDrawStyle::TCanvasApplyCssStyle(const char *styleName, TCanvas *cCanvas, Int_t verbose) {

  Bool_t status = false;
  Bool_t status2 = false;
  TString elementID = "";
  TString classID = "";
  TString objectID = "";
  TString property = "";
  TString localStyle = "";
  Int_t valueI;
  Float_t valueF;

  AliDrawStyle::GetIds(cCanvas, elementID, classID, objectID, localStyle, verbose);

  Int_t cWidth = 0;
  Int_t cHeight = 0;
  property = AliDrawStyle::GetValue(styleName, "width", elementID, classID, objectID, localStyle);
  if (property != "") cWidth = (Int_t) AliDrawStyle::GetNamedTypeAt(property, status, 0);
  property = AliDrawStyle::GetValue(styleName, "height", elementID, classID, objectID, localStyle);
  if (property != "") cHeight = (Int_t) AliDrawStyle::GetNamedTypeAt(property, status2, 0);
  if (cWidth * cHeight > 0 && status && status2) cCanvas->SetWindowSize((UInt_t) cWidth, (UInt_t) cHeight);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("fill-color"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetFillColor((Color_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("highlight-color"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetHighLightColor((Color_t) valueI);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-bottom"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetBottomMargin(valueF);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-top"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetTopMargin(valueF);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-left"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetLeftMargin(valueF);

  valueF = AliDrawStyle::PrepareValue(styleName, TString("margin-right"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetRightMargin(valueF);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("border-size"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetBorderSize((Short_t) valueI);

  valueI = (Int_t) AliDrawStyle::PrepareValue(styleName, TString("border-mode"), elementID, classID, objectID, localStyle, status, AliDrawStyle::GetPadNumber(), verbose);
  if (status) cCanvas->SetBorderMode((Short_t) valueI);

}

///
/// \brief Method get IDs from the cObject->GetName() and parse it into input reference.
///        Names of objects should according to css notation.
///        #### Example of name:
///        \code
///           TH1F *hisArray[nHis];
///           for (Int_t i=0; i<nHis; i++) {
///             hisArray[i] = new TH1F(TString::Format("his[%d].class(Raw,Error)",i).Data(),
///                                    TString::Format("his[%d].class(Raw,Error)",i).Data(), 100, -2, 2);
///           }
///        \endcode
///         Schema of parsing:
///
/// \param cObject    - input TObject
/// \param elementID  - cObject->GetName()           (in example - TH1)
/// \param classID    - values from .class(*)        (in example - Raw,Error)
/// \param objectID   - value before number or class (in example - his )
  void
  AliDrawStyle::GetIds(TObject *cObject, TString &elementID, TString &classID, TString &objectID, TString &localStyle, Int_t verbose) {

  elementID = cObject->ClassName();
  objectID = TString(cObject->GetName());
  TPRegexp classPat(".class[(].*?[)]");
  TPRegexp stylePat(".style[(].*?[)]");
  TString objPatStr = "";
  if ((objectID.Index(".class") < objectID.Index(".style")  || objectID.Index(".style") == -1) && objectID.Index(".class") > 0)
    objPatStr = "^.*[[]|^.*?[.]class|.*";
  else if ((objectID.Index(".class") > objectID.Index(".style") || objectID.Index(".class") == -1) && objectID.Index(".style") > 0)
    objPatStr = "^.*[[]|^.*?[.]style|.*";
  else
    objPatStr = "^.*[[]|.*";
  TPRegexp objPat(objPatStr);

  classID = objectID(classPat);
  localStyle = objectID(stylePat);
  Int_t classLength = classID.Index(")") - classID.Index("(") - 1;
  Int_t styleLength = localStyle.Index(")") - localStyle.Index("(") - 1;
  classID = classID(classID.Index("(") + 1, classLength);
  localStyle = localStyle(localStyle.Index("(") + 1, styleLength);
  objectID = objectID(objPat);
  objectID.ReplaceAll(".class", "");
  objectID.ReplaceAll(".style", "");
  if (objectID(objectID.Length() - 1) == '[')
    objectID = objectID(0,objectID.Length() - 1);
    if (verbose == 4)
      ::Info("AliDrawStyle", "Object with name \"%s\" was parsed via AliDrawStyle::GetIds() to elementID = \"%s\", classID = \"%s\", objectID = \"%s\", localStyle = \"%s\"", cObject->GetName(), elementID.Data(), classID.Data(), objectID.Data(), localStyle.Data());
}

/// \brief Applies style from css to all objects from Pad or Canvas.
///        In case if pad inherited from TCanvas will work recursively for all pads from input canvas.
/// \param pad       - Input TPad object. You can specify TCanvas in this case style will be apply recursively to all objects (TH1, TF1, TGraph) on pad.
/// \param styleName - Name of style specify in AliDrawStyle::RegisterCssStyle()
/*!
  ####  Example usage see in the beginnings of this file:
*/
void AliDrawStyle::ApplyCssStyle(TPad *pad, const char *styleName, Int_t verbose) {
  if (pad == nullptr) {
    ::Error("AliDrawStyle::ApplyCssStyle", "Pad doesn't exist");
    return;
  }

  TObject *cObj = nullptr;
  TList *oList = nullptr;
  TString elementID = "";
  TString objectID = "";
  TString classID = "";
  TString localStyle = "";
  TPRegexp numPat0("[[].*[]]");
  TPRegexp numPat1("[0-9]+");
  TString padNumStr;

  oList = pad->GetListOfPrimitives();
  GetIds(pad, elementID, classID, objectID, localStyle, verbose);

  if (elementID == "TCanvas") {
    AliDrawStyle::TCanvasApplyCssStyle(styleName, (TCanvas *) pad, verbose);
    pad->Modified();
    for (Int_t c = 0; c < oList->GetEntries() && oList->At(c)->InheritsFrom("TPad"); c++) {
      AliDrawStyle::ApplyCssStyle((TPad *) oList->At(c), styleName, verbose);
    }
  }

//fixme: for pads we can use indexes only from name
  if (elementID == "TPad") {
    padNumStr = TString(TString(pad->GetName())(numPat0))(numPat1);
    if (padNumStr != TString(""))
      AliDrawStyle::SetPadNumber(padNumStr.Atoi());
      AliDrawStyle::TPadApplyStyle(styleName, pad, verbose);
      AliDrawStyle::SetPadNumber(AliDrawStyle::GetPadNumber() + 1);
  }
  Int_t objNum = -1, hisCnt = 0, funcCnt = 0, graphCnt = 0, legendCnt = 0;
  for (Int_t k = 0; k < oList->GetEntries(); k++) {
    cObj = oList->At(k);
    if (TString(TString(cObj->GetName())(numPat0))(numPat1) != TString(""))
      objNum = TString(TString(TString(cObj->GetName())(numPat0))(numPat1)).Atoi();
    if (cObj->InheritsFrom("TH1")) {
      if (objNum >= 0) hisCnt = objNum;
      AliDrawStyle::TH1ApplyStyle(styleName, (TH1 *) cObj, hisCnt, verbose);
      hisCnt++;
    }
    if (cObj->InheritsFrom("TGraph")) {
      if (objNum >= 0) graphCnt = objNum;
      AliDrawStyle::TGraphApplyStyle(styleName, (TGraph *) cObj, graphCnt, verbose);
      graphCnt++;
    }
    if (cObj->InheritsFrom("TF1")) {
      if (objNum >= 0) funcCnt = objNum;
      AliDrawStyle::TF1ApplyStyle(styleName, (TF1 *) cObj, funcCnt, verbose);
      funcCnt++;
    }
    if (cObj->InheritsFrom("TLegend")) {
      if (objNum >= 0) legendCnt = objNum;
      AliDrawStyle::TLegendApplyStyle(styleName, (TLegend *) cObj, legendCnt, verbose);
      legendCnt++;
    }
  }
  pad->Modified();
}
