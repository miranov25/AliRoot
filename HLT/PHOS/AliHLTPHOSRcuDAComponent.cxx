/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliPHOSRcuDA1.h"
#include "AliHLTPHOSSharedMemoryInterface.h"
#include "AliHLTPHOSRcuDAComponent.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "TObjArray.h"

//#include <iostream>

/** @file   AliHLTPHOSRcuDAComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  A module calibration component for PHOS HLT, using the PHOS DA's
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

using namespace PhosHLTConst;

AliHLTPHOSRcuDAComponent gAliHLTPHOSRcuDAComponent;

AliHLTPHOSRcuDAComponent::AliHLTPHOSRcuDAComponent() : AliHLTPHOSRcuProperties(),
						       AliHLTCalibrationProcessor(),
						       fPhosEventCount(0),
						       fPHOSDAPtr(0),
						       fShmPtr(0) 
						       //    fTest(-2)
{
  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
}



AliHLTPHOSRcuDAComponent::~AliHLTPHOSRcuDAComponent() 
{
  if(fShmPtr)
    {
      delete fShmPtr;
      fShmPtr = 0;
    }
  if(fPHOSDAPtr)
    {
      delete fPHOSDAPtr;
      fPHOSDAPtr = 0;
    }
}


void AliHLTPHOSRcuDAComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCellEnergyDataType);
}


AliHLTComponentDataType AliHLTPHOSRcuDAComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkEmcCalibDataType;
}
  
                                   
void AliHLTPHOSRcuDAComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 0;
  inputMultiplier = 2;
}


AliHLTComponent* 
AliHLTPHOSRcuDAComponent::Spawn()
{
  return new AliHLTPHOSRcuDAComponent();
}

const char* 
AliHLTPHOSRcuDAComponent::GetComponentID()
{
  return "PhosRcuDAProcessor";  
}


Int_t
AliHLTPHOSRcuDAComponent::ScanArgument( Int_t argc, const char** argv)
{
  ScanArguments(argc, argv);
  return 0;
}


Int_t AliHLTPHOSRcuDAComponent::InitCalibration()
{  
  //CRAP PT just to get something working by 5 May 2008
  const int tmpModule = 2;
  fPHOSDAPtr = new AliPHOSRcuDA1(tmpModule ,GetRCUID()); 
  return 0;
}


Int_t AliHLTPHOSRcuDAComponent::DeinitCalibration()
{
  AliHLTComponentEventData dummyEvtData;
  AliHLTComponentTriggerData dummyTrgData;
  ShipDataToFXS(dummyEvtData, dummyTrgData); 


  if(fShmPtr)
    {
      delete fShmPtr;
      fShmPtr = 0;
    }
  if(fPHOSDAPtr)
    {
      delete fPHOSDAPtr;
      fPHOSDAPtr = 0;
    }
  return 0;
}



Int_t AliHLTPHOSRcuDAComponent::ProcessCalibration(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  fPhosEventCount ++;
  const  AliHLTComponentEventData eDta  = evtData;
  AliHLTComponentTriggerData  tDta =  trigData;

  UInt_t specification = 0;
  const AliHLTComponentBlockData* iter = 0;
  iter = GetFirstInputBlock( AliHLTPHOSDefinitions::fgkCellEnergyDataType | kAliHLTDataOriginPHOS);

  AliHLTPHOSRcuCellEnergyDataStruct* cellDataPtr = 0;
  Int_t xOffset = 0;
  Int_t zOffset = 0;
  Int_t module = -1;

  Float_t energyArray[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];
  Float_t timeArray[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];
  ResetArrays(energyArray, timeArray);

  while(iter != 0)
    {
      specification = specification|iter->fSpecification;
      cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      module = cellDataPtr->fModuleID;
      xOffset = cellDataPtr->fRcuX*N_XCOLUMNS_RCU;
      zOffset = cellDataPtr->fRcuZ*N_ZROWS_RCU;

      for(Int_t x = 0; x < N_XCOLUMNS_RCU; x++)
	{
	  for(Int_t z = 0; z < N_ZROWS_RCU; z++)
	    {
	      for(Int_t gain = 0; gain < N_GAINS; gain++)
		{
		  energyArray[x+xOffset][z+zOffset][gain] = cellDataPtr->fValidData[x][z][gain].fEnergy;
		  timeArray[x+xOffset][z+zOffset][gain] = cellDataPtr->fValidData[x][z][gain].fTime;
		}
	    }
	}
      iter = GetNextInputBlock(); 
    }
  
  fPHOSDAPtr->FillHistograms(energyArray, timeArray);

  ResetArrays(energyArray, timeArray);

  return 0; 
}

  
Int_t 
AliHLTPHOSRcuDAComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) 
{
  Char_t filename[200];

  for(int i=0; i < 200; i++)
    {
      filename[i] = 0;
    }
  const TObjArray *calibPtr = fPHOSDAPtr->GetHistoContainer();
  sprintf(filename, "/home/perthi/hlt/rundir/test/outdata/%s.root", fPHOSDAPtr->GetName() );
  TFile *outFile =  new TFile(filename, "recreate");
  calibPtr->Write(); 
  outFile->Close();
  PushToFXS( (TObject*)fPHOSDAPtr->GetHistoContainer(), "PHOS",  filename);
  cout << "Finnished pushing data to HLT FXS" << endl;
  return 0;
}  


void
AliHLTPHOSRcuDAComponent::ResetArrays(Float_t e[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS], Float_t t[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS])
{
  for(Int_t x = 0; x < N_XCOLUMNS_RCU; x++)
    {
      for(Int_t z = 0; z < N_ZROWS_RCU; z++)
	{
	  for(Int_t gain = 0; gain < N_GAINS; gain++)
	    {
	      e[x][z][gain] = 0;
	      t[x][z][gain] = 0;
	    }
	}
    }
}
