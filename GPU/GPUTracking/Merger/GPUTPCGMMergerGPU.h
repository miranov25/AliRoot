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

/// \file GPUTPCGMMergerGPU.h
/// \author David Rohr

#ifndef GPUTPCGMMERGERGPUCA_H
#define GPUTPCGMMERGERGPUCA_H

#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCGMMergerTrackFit : public GPUKernelTemplate
{
 public:
  GPUhdi() CONSTEXPR static GPUDataTypes::RecoStep GetRecoStep() { return GPUDataTypes::RecoStep::TPCMerging; }
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  typedef GPUTPCGMMerger processorType;
  GPUhdi() static processorType* Processor(GPUConstantMem& processors)
  {
    return &processors.tpcMerger;
  }
  template <int iKernel = 0>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUTPCSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerUnpack
{
};
class GPUTPCGMMergerMergeWithin
{
};
class GPUTPCGMMergerMergeSlices
{
};
class GPUTPCGMMergerMergeCE
{
};
class GPUTPCGMMergerCollect
{
};
class GPUTPCGMMergerClusters
{
};
class GPUTPCGMMergerCopyToGPU
{
};
class GPUTPCGMMergerCopyToHost
{
};
class GPUTPCGMMergerFinalize
{
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
