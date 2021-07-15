#ifndef ALITPCPIDRESPONSE_H
#define ALITPCPIDRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TPC PID class
// A very naive design... Should be made better by the detector experts...
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
// With many additions and modifications suggested by
//      Alexander Kalweit, GSI, alexander.philipp.kalweit@cern.ch
//      Dariusz Miskowiec, GSI, D.Miskowiec@gsi.de
// ...and some modifications by
//      Mikolaj Krzewicki, GSI, mikolaj.krzewicki@cern.ch
// ...and some modifications plus eta correction functions by
//      Benjamin Hess, University of Tuebingen, bhess@cern.ch
//-------------------------------------------------------
#include <Rtypes.h>

#include <TNamed.h>
#include <TVectorF.h>
#include <TObjArray.h>
#include <TF1.h>
#include <TString.h>

#include "AliPID.h"
#include "AliVTrack.h"

class TH2D;
class TSpline3;
class AliOADBContainer;
class AliNDLocalRegression;

class AliTPCPIDResponse: public TNamed {
public:
  AliTPCPIDResponse();
  //TODO Remove? AliTPCPIDResponse(const Double_t *param);
  AliTPCPIDResponse(const AliTPCPIDResponse&);
  AliTPCPIDResponse& operator=(const AliTPCPIDResponse&);
  virtual ~AliTPCPIDResponse();

  enum EChamberStatus {
    kChamberOff=0,
    kChamberHighGain=1,
    kChamberLowGain=2,
    kChamberInvalid=3
  };
  
  enum ETPCgainScenario {
    kDefault= 0,
    kALLhigh = 1,
    kOROChigh = 2,
    kGainScenarioInvalid = 3
  };

  static const Int_t fgkNumberOfParticleSpecies=AliPID::kSPECIESC;
  static const Int_t fgkNumberOfGainScenarios=3;
  static const Int_t fgkNumberOfdEdxSourceScenarios=3;

  enum ETPCdEdxSource {
    kdEdxDefault=0,        // use combined dEdx from IROC+OROC (assumes ideal detector)
    kdEdxOROC=1,       // use only OROC
    kdEdxHybrid=2,   // Use IROC+OROC dEdx only where IROCS are good (high gain), otherwise fall back to OROC only
    kdEdxInvalid=3     //invalid
  };

  enum ETPCdEdxType {
    kdEdxTrack=0,
    kdEdxInfo=1
  };

  enum ETPCPileupCorrectionStrategy {
    kPileupCorrectionInExpectedSignal = 0, // modifiy the expected signal as done for the eta and multiplicity correction
    kPileupCorrectionInTrackSignal = 1, // directly modify the dEdx
    kNoPileupCorrection = 2 // don't do pileup correction, this takes preceedence over the flag in GetExpectedSigma, etc
  };

  enum EMultiplicityEstimator {
    kNumberOfESDTracks = 0, // get multiplicty from AliVEvent::GetNumberOfESDTracks()
    kNTPCTrackBeforeClean = 1 // get multiplicty from AliVEvent::GetNTPCTrackBeforeClean()
  };

  void SetSigma(Float_t res0, Float_t resN2);
  void SetBetheBlochParameters(Double_t kp1,
                               Double_t kp2,
                               Double_t kp3,
                               Double_t kp4,
                               Double_t kp5
                               );
  //Better prevent user from setting fMIP != 50. because fMIP set fix to 50 for much other code:
  void SetMip(Float_t mip) { fMIP = mip; } // Set overall normalisation; mean dE/dx for MIP
  Double_t Bethe(Double_t bg) const;
  void SetUseDatabase(Bool_t useDatabase) { fUseDatabase = useDatabase;}
  Bool_t GetUseDatabase() const { return fUseDatabase;}
  
  void SetResponseFunction(AliPID::EParticleType type, TObject * const o) { fResponseFunctions.AddAt(o,(Int_t)type); }
  const TObject * GetResponseFunction(AliPID::EParticleType type) const { return fResponseFunctions.At((Int_t)type); }
  void SetVoltage(Int_t n, Float_t v) {fVoltageMap[n]=v;}
  void SetVoltageMap(const TVectorF& a) {fVoltageMap=a;} //resets ownership, ~ will not delete contents
  Float_t GetVoltage(Int_t n) const {return fVoltageMap[n];}
  void SetLowGainIROCthreshold(Float_t v) {fLowGainIROCthreshold=v;}
  void SetBadIROCthreshold(Float_t v) {fBadIROCthreshhold=v;}
  void SetLowGainOROCthreshold(Float_t v) {fLowGainOROCthreshold=v;}
  void SetBadOROCthreshold(Float_t v) {fBadOROCthreshhold=v;}
  void SetMaxBadLengthFraction(Float_t f) {fMaxBadLengthFraction=f;}

  void SetMagField(Double_t mf) { fMagField=mf; }
  
  const TH2D* GetEtaCorrMap() const { return fhEtaCorr; };
  Bool_t SetEtaCorrMap(TH2D* hMap);
  
  Double_t GetTrackTanTheta(const AliVTrack *track) const;
  
  const TH2D* GetSigmaPar1Map() const { return fhEtaSigmaPar1; };
  Double_t GetSigmaPar0() const { return fSigmaPar0; };
  Bool_t SetSigmaParams(TH2D* hSigmaPar1Map, Double_t sigmaPar0);
  
  Double_t GetSigmaPar1(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource = kdEdxDefault) const;
  
  Bool_t GetIsNewPbPb() const {return fIsNewPbPbParam; };
 
  const TF1* GetMultSlopeFunction() const  { return fCorrFuncSlope; };
  void SetParameterMultSlope(Int_t parIndex, Double_t parValue)  
      { if (fCorrFuncSlope) fCorrFuncSlope->SetParameter(parIndex, parValue); };
      
  const TF1* GetMultCurvFunction() const  { return fCorrFuncCurv; };
  void SetParameterMultCurv(Int_t parIndex, Double_t parValue)  
      { if (fCorrFuncCurv) fCorrFuncCurv->SetParameter(parIndex, parValue); };      
  
  const TF1* GetMultiplicityCorrectionFunction() const  { return fCorrFuncMultiplicity; };
  void SetParameterMultiplicityCorrection(Int_t parIndex, Double_t parValue)  
      { if (fCorrFuncMultiplicity) fCorrFuncMultiplicity->SetParameter(parIndex, parValue); };
  
  const TF1* GetMultiplicityCorrectionFunctionTanTheta() const  { return fCorrFuncMultiplicityTanTheta; };
  void SetParameterMultiplicityCorrectionTanTheta(Int_t parIndex, Double_t parValue)  
      { if (fCorrFuncMultiplicityTanTheta) fCorrFuncMultiplicityTanTheta->SetParameter(parIndex, parValue); };

  const TF1* GetMultiplicitySigmaCorrectionFunction() const  { return fCorrFuncSigmaMultiplicity; };
  void SetParameterMultiplicitySigmaCorrection(Int_t parIndex, Double_t parValue)  
      { if (fCorrFuncSigmaMultiplicity) fCorrFuncSigmaMultiplicity->SetParameter(parIndex, parValue); };
  
  void ResetMultiplicityCorrectionFunctions(); 
  
  void SetCurrentEventMultiplicity(Int_t value) { fCurrentEventMultiplicity = value;  };
  Int_t GetCurrentEventMultiplicity() const { return fCurrentEventMultiplicity; };

  EMultiplicityEstimator GetMultiplicityEstimator() const { return fMultiplityEstimator; }
  void SetMultiplicityEstimator(EMultiplicityEstimator estimator) { fMultiplityEstimator = estimator; }

  Double_t GetEtaCorrection(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource = kdEdxDefault) const;
  Double_t GetMultiplicityCorrection(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource = kdEdxDefault) const;
  Double_t GetMultiplicitySigmaCorrection(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource = kdEdxDefault) const;

  Double_t GetMultiplicityCorrectedTrackdEdx(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource = kdEdxDefault) const;
  Double_t GetEtaAndMultiplicityCorrectedTrackdEdx(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource = kdEdxDefault) const;
  Double_t GetEtaCorrectedTrackdEdx(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource = kdEdxDefault) const;
  Double_t GetCorrectedTrackdEdx(const AliVTrack *track, AliPID::EParticleType species, Bool_t applyEtaCorrection, Bool_t applyMultiplicityCorrection, Bool_t applyPileupCorrection, ETPCdEdxSource dedxSource = kdEdxDefault);
  
  // Fast functions for expert use only
  Double_t GetEtaCorrectionFast(const AliVTrack *track, Double_t dEdxSplines) const;
  
  Double_t GetMultiplicityCorrectionFast(const AliVTrack *track, Double_t dEdxExpected, Int_t multiplicity) const;
  
  Double_t GetMultiplicitySigmaCorrectionFast(Double_t dEdxExpected, Int_t multiplicity) const;
  
  Double_t GetSigmaPar1Fast(const AliVTrack *track, AliPID::EParticleType species,
                            Double_t dEdx, const TSpline3* responseFunction) const;
  
  //NEW
  void SetSigma(Float_t res0, Float_t resN2, ETPCgainScenario gainScenario );
  Double_t GetExpectedSignal( const AliVTrack* track,
                              AliPID::EParticleType species,
                              ETPCdEdxSource dedxSource = kdEdxDefault,
                              Bool_t correctEta = kFALSE,
                              Bool_t correctMultiplicity = kFALSE,
                              Bool_t usePileupCorrection = kFALSE,
                              Bool_t useQPtTglCorrection = kFALSE,
                              Bool_t useEnergyLossCorrection = kFALSE) const;
  Double_t GetExpectedSigma( const AliVTrack* track, 
                             AliPID::EParticleType species,
                             ETPCdEdxSource dedxSource = kdEdxDefault,
                             Bool_t correctEta = kFALSE,
                             Bool_t correctMultiplicity = kFALSE,
                             Bool_t usePileupCorrection = kFALSE,
                              Bool_t useQPtTglCorrection = kFALSE,
                              Bool_t useEnergyLossCorrection = kFALSE) const;
  Float_t GetNumberOfSigmas( const AliVTrack* track,
                             AliPID::EParticleType species,
                             ETPCdEdxSource dedxSource = kdEdxDefault,
                             Bool_t correctEta = kFALSE,
                             Bool_t correctMultiplicity = kFALSE,
                             Bool_t usePileupCorrection = kFALSE,
                            Bool_t useQPtTglCorrection = kFALSE,
                              Bool_t useEnergyLossCorrection = kFALSE) const;
  
  Float_t GetSignalDelta( const AliVTrack* track,
                          AliPID::EParticleType species,
                          ETPCdEdxSource dedxSource = kdEdxDefault,
                          Bool_t correctEta = kFALSE,
                          Bool_t correctMultiplicity = kFALSE,
                          Bool_t usePileupCorrection = kFALSE,
                          Bool_t useQPtTglCorrection = kFALSE,
                          Bool_t useEnergyLossCorrection = kFALSE,
                          Bool_t ratio = kFALSE) const;
  
  void SetResponseFunction(TObject* o,
                           AliPID::EParticleType type,
                           ETPCgainScenario gainScenario);
  void Print(Option_t* option="") const;
  TSpline3* GetResponseFunction( AliPID::EParticleType species,
                                 ETPCgainScenario gainScenario ) const;
  TSpline3* GetResponseFunction( const AliVTrack* track,
                                 AliPID::EParticleType species,
                                 ETPCdEdxSource dedxSource = kdEdxDefault) const;
  Bool_t ResponseFunctiondEdxN(const AliVTrack* track, 
                               AliPID::EParticleType species,
                               ETPCdEdxSource dedxSource,
                               Double_t& dEdx, Int_t& nPoints, ETPCgainScenario& gainScenario, TSpline3** responseFunction) const;
  Bool_t sectorNumbersInOut(Double_t* trackPositionInner,
                            Double_t* trackPositionOuter,
                            Float_t& phiIn, Float_t& phiOut, 
                            Int_t& in, Int_t& out ) const;
  AliTPCPIDResponse::EChamberStatus TrackStatus(const AliVTrack* track, Int_t layer) const;
  Float_t MaxClusterRadius(const AliVTrack* track) const;
  Bool_t TrackApex(const AliVTrack* track, Float_t magField, Double_t position[3]) const;
  static const char* GainScenarioName(Int_t n) {return fgkGainScenarioName[(n>fgkNumberOfGainScenarios)?fgkNumberOfGainScenarios:n];}
  Int_t ResponseFunctionIndex( AliPID::EParticleType species,
                               ETPCgainScenario gainScenario ) const;
  void ResetSplines();

  //OLD
  Double_t GetExpectedSignal(Float_t mom,
                     AliPID::EParticleType n=AliPID::kKaon) const;
  Double_t GetExpectedSigma(Float_t mom, Int_t nPoints,
                            AliPID::EParticleType n=AliPID::kKaon) const;
  Float_t  GetNumberOfSigmas(Float_t mom, 
                             Float_t dEdx, 
			     Int_t nPoints,
                             AliPID::EParticleType n=AliPID::kKaon) const {
    //
    // Deprecated function (for backward compatibility). Please use 
    // GetNumberOfSigmas(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource,
    // Bool_t correctEta, Bool_t correctMultiplicity)
    // instead!TODO
    //
    
    Double_t bethe=GetExpectedSignal(mom,n);
    Double_t sigma=GetExpectedSigma(mom,nPoints,n);
    return (dEdx-bethe)/sigma;
  }

  Double_t GetMIP() const { return fMIP;} 
  Float_t  GetRes0()  const { return fRes0[0];  }
  Float_t  GetResN2() const { return fResN2[0]; }
  Float_t  GetRes0(ETPCgainScenario s)  const { return fRes0[s];  }
  Float_t  GetResN2(ETPCgainScenario s) const { return fResN2[s]; }

  Bool_t   RegisterSpline(const char * name, Int_t index);
  Double_t EvaldEdxSpline(Double_t bg,Int_t entry);
  static   Double_t SEvaldEdx(Double_t bg,Int_t entry){ return (fgInstance!=0)? fgInstance->EvaldEdxSpline(bg,entry):0;};

  //===| Pileup correction |====================================================
  void SetEventPileupProperties(Double_t shift, Double_t pileup, Double_t mult) { fEventPileupProperties[0] = shift; fEventPileupProperties[1] = pileup; fEventPileupProperties[2] = mult; }
  Float_t GetPileUpProperties(UInt_t i ) {return fEventPileupProperties[i%3];}
  void SetPileupCorrectionStrategy(ETPCPileupCorrectionStrategy strategy) { fPileupCorrectionStrategy = strategy; }
  // pileup correction
  void SetPileupCorrectionObject(AliNDLocalRegression* correction) { fPileupCorrection = correction; }
  const AliNDLocalRegression* GetPileupCorrectionObject() const { return fPileupCorrection; }
  Bool_t IsPileupCorrectionRequested() const { return fPileupCorrectionRequested; }
  Double_t GetPileupCorrectionValue(const AliVTrack* track) const;
  static AliNDLocalRegression* GetPileupCorrectionFromFile(const TString fileName);
  //qPtTgl correction
  void SetQPtTglCorrectionObject(AliNDLocalRegression* correction) { fQPtTglCorrection = correction; }
  const AliNDLocalRegression* GetQPtTglCorrectionObject() const { return fQPtTglCorrection; }
  Bool_t IsQPtTglCorrectionRequested() const { return fQPtTglCorrectionRequested; }
  Double_t GetQPtTglCorrectionValue(const AliVTrack* track) const;
  static AliNDLocalRegression* GetQPtTglCorrectionFromFile(const TString fileName);
  // energy loss correction
  void SetEnergyLossCorrection(Int_t gasType, Float_t norm){fGasType=gasType; fEnergyLossFactor=norm;}
  Bool_t IsEnergyLossCorrectionRequested() const { return fEnergyLossCorrectionRequested; }
  Double_t GetEnergyLossCorrectionValue(const AliVTrack* track, Int_t pidCode) const;
  //
  //===| dEdx type functions |==================================================
  void SetdEdxType(ETPCdEdxType dEdxType, Int_t dEdxChargeType=0, Int_t dEdxWeightType=0, Double_t dEdxIROCweight=1., Double_t dEdxOROCmedWeight=1., Double_t dEdxOROClongWeight=1.) {
    fdEdxType=dEdxType; fdEdxChargeType=dEdxChargeType; fdEdxWeightType=dEdxWeightType; fIROCweight=dEdxIROCweight; fOROCmedWeight=dEdxOROCmedWeight; fOROClongWeight=dEdxOROClongWeight; }
  ETPCdEdxType      GetdEdxType()        const { return fdEdxType;       }
  Int_t             GetdEdxChargeType()  const { return fdEdxChargeType; }
  Int_t             GetdEdxWeightType()  const { return fdEdxWeightType; }
  Double_t          GetIROCweight()      const { return fIROCweight;     }
  Double_t          GetOROCmedWeight()   const { return fOROCmedWeight;  }
  Double_t          GetOROClongWeight()  const { return fOROClongWeight; }

  Double_t GetTrackdEdx(const AliVTrack* track) const;

  //===| New resolution parametrization |=======================================
  void SetResolutionParametrization(TF1* fun) { fSigmaParametrization = fun; }
  TF1* GetResolutionParametrization() const { return fSigmaParametrization; }

  void SetMultiplicityNormalization(Double_t norm) { fMultiplicityNormalization = norm; }
  Double_t GetMultiplicityNormalization() const { return fMultiplicityNormalization; }

  /// index PID estimator as in AliTPCdEdxInfo: 0-IROC, 1-OROCmedium, 2-OROClong, 3-OROCall, 4-FullTrack
  Double_t GetExpectedSigmaTF1(const AliVTrack* track,
                               AliPID::EParticleType species,
                               Int_t dEdxType = 4) const;

  void GetTF1ParametrizationValues(Double_t values[7],
                                   const AliVTrack* track,
                                   AliPID::EParticleType species,
                                   Int_t dEdxType = 4) const;

  //===| Initialisation |=======================================================
  Bool_t InitFromOADB(const Int_t run, Int_t pass, TString passName,
                      const char* oadbFile="$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCPIDResponseOADB.root",
                      Bool_t initMultiplicityCorrection=kTRUE);

  Bool_t SetSplinesFromArray                (const TObjArray* arrSplines);
  Bool_t SetMultiplicityCorrectionFromString(const TString& multiplicityData);
  Bool_t SetdEdxTypeFromString              (const TString& dEdxTypeSet);
  Bool_t SetdEdxResolutionFromString        (const TString& dEdxTypeSet);

  const TString& GetRecoPassNameUsed() { return  fRecoPassNameUsed; }

  //===| Helpers |==============================================================
  static TString GetChecksum(const TObject* obj);
  static TObjArray* GetMultiplicityCorrectionArrayFromString(const TString& corrections);
  //
  static Double_t BetheBlochAleph(Double_t bg, Double_t kp1, Double_t kp2, Double_t kp3, Double_t kp4, Double_t kp5);
  static double   sigmadEdxPt(double p, double PID, double resol=0.01 );
  static double GetPtOutHelix(Float_t ptIn, Double_t tgl, Float_t mass, Double_t rIn=84, Double_t rOut=245,  Int_t type=0, Float_t bz=5,Float_t fraction=1);
  static double dEdxRationR0R1Helix(Float_t ptIn, Float_t tgl,  Float_t mass, Double_t rIn0=93,Float_t rOut0=245, Double_t rIn1=93,Float_t rOut1=245,  Int_t type=0, Float_t bz=5, Int_t nStep=20);
  static double dEdxMeanToInHelix(Float_t ptIn, Float_t tgl,  Float_t mass, Double_t rIn=93,Float_t rOut=245,  Int_t type=0, Float_t bz=5, Int_t nStep=10, Float_t scale=1, Bool_t useInnerPt=kFALSE);

protected:
  Double_t GetExpectedSignal(const AliVTrack* track,
                             AliPID::EParticleType species,
                             Double_t dEdx,
                             const TSpline3* responseFunction,
                             Bool_t correctEta,
                             Bool_t correctMultiplicity,
                             Bool_t usePileupCorrection,
                             Bool_t useQPtTglCorrection = kFALSE,
                              Bool_t useEnergyLossCorrection = kFALSE) const;
  
  Double_t GetExpectedSigma(const AliVTrack* track, 
                            AliPID::EParticleType species,
                            ETPCgainScenario gainScenario,
                            Double_t dEdx,
                            Int_t nPoints,
                            const TSpline3* responseFunction,
                            Bool_t correctEta,
                            Bool_t correctMultiplicity,
                            Bool_t usePileupCorrection,
                            Bool_t useQPtTglCorrection = kFALSE,
                              Bool_t useEnergyLossCorrection = kFALSE) const;
  //
  // function for numberical debugging 0 registed splines can be used in the TFormula and tree visualizations
  //
private:
  Float_t fMIP;          // dEdx for MIP
  Float_t fRes0[fgkNumberOfGainScenarios];  // relative dEdx resolution  rel sigma = fRes0*sqrt(1+fResN2/npoint)
  Float_t fResN2[fgkNumberOfGainScenarios]; // relative Npoint dependence rel  sigma = fRes0*sqrt(1+fResN2/npoint)

  Double_t fKp1;   // Parameters
  Double_t fKp2;   //    of
  Double_t fKp3;   // the ALEPH
  Double_t fKp4;   // Bethe-Bloch
  Double_t fKp5;   // formula

  Bool_t fUseDatabase; // flag if fine-tuned database-response or simple ALEPH BB should be used
  
  TObjArray fResponseFunctions; //! ObjArray of response functions individually for each particle
  AliOADBContainer* fOADBContainer; //! OADB container with response functions
  AliNDLocalRegression* fPileupCorrection; // pileup correction object
  AliNDLocalRegression* fQPtTglCorrection; // pileup correction object
  Bool_t fEnergyLossCorrectionRequested;   //  flag fEnergyLossCorrectionRequested
  Int_t   fGasType;                        // type of the gas (0-Argon, 1-Neon)
  Float_t  fEnergyLossFactor;              // energy loss facto - fit parameters in respect to tabulated value
  TVectorF fVoltageMap; //!stores a map of voltages wrt nominal for all chambers
  Float_t fLowGainIROCthreshold;  //voltage threshold below which the IROC is considered low gain
  Float_t fBadIROCthreshhold;     //voltage threshold for bad IROCS
  Float_t fLowGainOROCthreshold;  //voltage threshold below which the OROC is considered low gain
  Float_t fBadOROCthreshhold;     //voltage threshold for bad OROCS
  Float_t fMaxBadLengthFraction;  //the maximum allowed fraction of track length in a bad sector.

  Int_t sectorNumber(Double_t phi) const;

  Double_t fMagField;  //! Magnetic field

  static const char* fgkGainScenarioName[fgkNumberOfGainScenarios+1];

  TH2D* fhEtaCorr; //! Map for TPC eta correction
  TH2D* fhEtaSigmaPar1; //! Map for parameter 1 of the IROCdEdx sigma parametrisation
  
  Double_t fSigmaPar0; // Parameter 0 of the dEdx sigma parametrisation
  
  EMultiplicityEstimator fMultiplityEstimator; // Multiplicity estimator to use
  Int_t fCurrentEventMultiplicity; // Multiplicity of the current event

  Double_t fEventPileupProperties[3]; //! current event properties: shift, pileup, multiplicity
  Bool_t fIsNewPbPbParam;
  TF1* fCorrFuncSlope;
  TF1* fCorrFuncCurv;
  TF1* fCorrFuncMultiplicity; //! Function to correct for the multiplicity dependence of the TPC dEdx
  TF1* fCorrFuncMultiplicityTanTheta; //! Function to correct the additional tanTheta dependence of the multiplicity dependence of the TPC dEdx
  TF1* fCorrFuncSigmaMultiplicity; //! Function to correct for the multiplicity dependence of the TPC dEdx resolution

  // dEdx type information
  ETPCdEdxType     fdEdxType;         // source for dEdx information to use
  Int_t            fdEdxChargeType;   // charge type to use for dEdx calculation from AliTPCdEdxInfo
  Int_t            fdEdxWeightType;   // tracklet weighting type to use for dEdx calculation from AliTPCdEdxInfo
  Double_t         fIROCweight;       // IROC weight to use for dEdx calculation from AliTPCdEdxInfo
  Double_t         fOROCmedWeight;    // OROC medium pad size weight to use for dEdx calculation from AliTPCdEdxInfo
  Double_t         fOROClongWeight;   // OROC long pad size weight to use for dEdx calculation from AliTPCdEdxInfo

  ETPCPileupCorrectionStrategy fPileupCorrectionStrategy; // Pileup correction strategy
  Bool_t fPileupCorrectionRequested; // If pileup correction was configured in the OADB object
  Bool_t fQPtTglCorrectionRequested; // If pileup correction was configured in the OADB object

  //===| New resolution parametrization |=======================================
  TF1*     fSigmaParametrization;      // Resolution parametrization
  Double_t fMultiplicityNormalization; // Value for the multiplicity normalisation in the sigma parametrization

  // Information on reconstruction data used
  TString fRecoPassNameUsed;          //! Name or number of the actually used reconstruction pass
  //
  //
  static AliTPCPIDResponse*   fgInstance;     //! Instance of this class (singleton implementation)
  TObjArray                   fSplineArray;   //array of registered splines
  ClassDef(AliTPCPIDResponse, 7)   // TPC PID class
};


#endif


