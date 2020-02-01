//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUSettings.cxx
/// \author David Rohr

#include "GPUSettings.h"
#include "GPUDef.h"
#include "GPUDataTypes.h"
#include <cstring>

using namespace GPUCA_NAMESPACE::gpu;

void GPUSettingsRec::SetDefaults()
{
  HitPickUpFactor = 2.;
  NeighboursSearchArea = 3.f;
  ClusterError2CorrectionY = 1.f;
  ClusterError2CorrectionZ = 1.f;
  MinNTrackClusters = -1;
  MaxTrackQPt = 1.f / GPUCA_MIN_TRACK_PT_DEFAULT;
  NWays = 1;
  NWaysOuter = false;
  RejectMode = 5;
  GlobalTracking = true;
  SearchWindowDZDR = 0.f;
  TrackReferenceX = 1000.f;
  NonConsecutiveIDs = false;
  ForceEarlyTPCTransform = -1;
  DisableRefitAttachment = 0;
  dEdxTruncLow = 2;
  dEdxTruncHigh = 77;
  tpcRejectionMode = GPUSettings::RejectionStrategyA;
  tpcRejectQPt = 1.f / 0.05f;
  tpcCompressionModes = GPUSettings::CompressionFull;
  tpcCompressionSortOrder = GPUSettings::SortPad;
  tpcSigBitsCharge = 4;
  tpcSigBitsWidth = 3;
  tpcZSthreshold = 2;
  fwdTPCDigitsAsClusters = 0;
}

void GPUSettingsEvent::SetDefaults()
{
  solenoidBz = -5.00668;
  constBz = 0;
  homemadeEvents = 0;
  continuousMaxTimeBin = 0;
  needsClusterer = 0;
}

void GPUSettingsProcessing::SetDefaults()
{
  deviceType = GPUDataTypes::DeviceType::CPU;
  forceDeviceType = true;
}

void GPUSettingsDeviceProcessing::SetDefaults()
{
  nThreads = 1;
  deviceNum = -1;
  platformNum = -1;
  globalInitMutex = false;
  gpuDeviceOnly = false;
  nDeviceHelperThreads = 2;
  debugLevel = -1;
  deviceTimers = false;
  debugMask = -1;
  comparableDebutOutput = true;
  resetTimers = 1;
  eventDisplay = nullptr;
  runQA = false;
  runCompressionStatistics = false;
  stuckProtection = 0;
  memoryAllocationStrategy = 0;
  keepAllMemory = false;
  nStreams = 8;
  trackletConstructorInPipeline = true;
  trackletSelectorInPipeline = false;
  forceMemoryPoolSize = 0;
  nTPCClustererLanes = 3;
}
