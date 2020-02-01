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

/// \file GPUTPCNeighboursFinder.h
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#ifndef GPUTPCNEIGHBOURSFINDER_H
#define GPUTPCNEIGHBOURSFINDER_H

#include "GPUTPCDef.h"
#include "GPUTPCRow.h"
#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
MEM_CLASS_PRE()
class GPUTPCTracker;

/**
 * @class GPUTPCNeighboursFinder
 *
 */
class GPUTPCNeighboursFinder : public GPUKernelTemplate
{
 public:
  MEM_CLASS_PRE()
  class GPUTPCSharedMemory
  {
    friend class GPUTPCNeighboursFinder;

   public:
#if !defined(GPUCA_GPUCODE)
    GPUTPCSharedMemory() : mNHits(0), mUpNHits(0), mDnNHits(0), mUpDx(0), mDnDx(0), mUpTx(0), mDnTx(0), mIRow(0), mIRowUp(0), mIRowDn(0), mRow(), mRowUp(), mRowDown()
    {
    }

    GPUTPCSharedMemory(const GPUTPCSharedMemory& /*dummy*/) : mNHits(0), mUpNHits(0), mDnNHits(0), mUpDx(0), mDnDx(0), mUpTx(0), mDnTx(0), mIRow(0), mIRowUp(0), mIRowDn(0), mRow(), mRowUp(), mRowDown() {}
    GPUTPCSharedMemory& operator=(const GPUTPCSharedMemory& /*dummy*/) { return *this; }
#endif //! GPUCA_GPUCODE

   protected:
    int mNHits;   // n hits
    int mUpNHits; // n hits in the next row
    int mDnNHits; // n hits in the prev row
    float mUpDx;  // x distance to the next row
    float mDnDx;  // x distance to the previous row
    float mUpTx;  // normalized x distance to the next row
    float mDnTx;  // normalized x distance to the previous row
    int mIRow;    // row number
    int mIRowUp;  // next row number
    int mIRowDn;  // previous row number
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0
    float2 mA[GPUCA_THREAD_COUNT_FINDER][GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP]; // temp memory
    calink mB[GPUCA_THREAD_COUNT_FINDER][GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP];
#endif
    MEM_LG(GPUTPCRow)
    mRow, mRowUp, mRowDown;
  };

  typedef GPUconstantref() MEM_GLOBAL(GPUTPCTracker) processorType;
  GPUhdi() CONSTEXPR static GPUDataTypes::RecoStep GetRecoStep() { return GPUCA_RECO_STEP::TPCSliceTracking; }
  MEM_TEMPLATE()
  GPUhdi() static processorType* Processor(MEM_TYPE(GPUConstantMem) & processors)
  {
    return processors.tpcTrackers;
  }
  template <int iKernel = 0>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUTPCSharedMemory) & smem, processorType& tracker);
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCNEIGHBOURSFINDER_H
