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

/// \file CfUtils.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_CF_UTILS_H
#define O2_GPU_CF_UTILS_H

#include "clusterFinderDefs.h"
#include "GPUCommonAlgorithm.h"
#include "Array2D.h"
#include "CfConsts.h"
#include "GPUTPCClusterFinderKernels.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class CfUtils
{

 public:
  static GPUdi() bool isAtEdge(const deprecated::Digit* d)
  {
    return (d->pad < 2 || d->pad >= TPC_PADS_PER_ROW - 2);
  }

  static GPUdi() bool innerAboveThreshold(uchar aboveThreshold, ushort outerIdx)
  {
    return aboveThreshold & (1 << CfConsts::OuterToInner[outerIdx]);
  }

  static GPUdi() bool innerAboveThresholdInv(uchar aboveThreshold, ushort outerIdx)
  {
    return aboveThreshold & (1 << CfConsts::OuterToInnerInv[outerIdx]);
  }

  static GPUdi() ushort partition(GPUTPCClusterFinderKernels::GPUTPCSharedMemory& smem, ushort ll, bool pred, ushort partSize, ushort* newPartSize)
  {
    bool participates = ll < partSize;

    ushort lpos = work_group_scan_inclusive_add((int)(!pred && participates));

    ushort part = work_group_broadcast(lpos, SCRATCH_PAD_WORK_GROUP_SIZE - 1);

    lpos -= 1;
    ushort pos = (participates && !pred) ? lpos : part;

    *newPartSize = part;
    return pos;
  }

// TODO: replace macros by templates
#if defined(GPUCA_GPUCODE)
#define DECL_FILL_SCRATCH_PAD(type, accessFunc)                                     \
  static GPUdi() void fillScratchPad_##type(                                        \
    GPUglobalref() const type* chargeMap,                                           \
    uint wgSize,                                                                    \
    uint elems,                                                                     \
    ushort ll,                                                                      \
    uint offset,                                                                    \
    uint N,                                                                         \
    GPUconstexprref() const Delta2* neighbors,                                      \
    GPUsharedref() const ChargePos* posBcast,                                       \
    GPUsharedref() type* buf)                                                       \
  {                                                                                 \
    GPUbarrier();                                                                   \
    ushort x = ll % N;                                                              \
    ushort y = ll / N;                                                              \
    Delta2 d = neighbors[x + offset];                                               \
    Delta dp = d.x;                                                                 \
    Delta dt = d.y;                                                                 \
    LOOP_UNROLL_ATTR for (unsigned int i = y; i < wgSize; i += (elems / N))         \
    {                                                                               \
      ChargePos readFrom = posBcast[i];                                             \
      uint writeTo = N * i + x;                                                     \
      buf[writeTo] = accessFunc(chargeMap, readFrom.gpad + dp, readFrom.time + dt); \
    }                                                                               \
    GPUbarrier();                                                                   \
  }

#define DECL_FILL_SCRATCH_PAD_COND(type, accessFunc, expandFunc, nameAppendix, null) \
  static GPUdi() void fillScratchPad##nameAppendix##_##type(                         \
    GPUglobalref() const type* chargeMap,                                            \
    uint wgSize,                                                                     \
    uint elems,                                                                      \
    ushort ll,                                                                       \
    uint offset,                                                                     \
    uint N,                                                                          \
    GPUconstexprref() const Delta2* neighbors,                                       \
    GPUsharedref() const ChargePos* posBcast,                                        \
    GPUsharedref() const uchar* aboveThreshold,                                      \
    GPUsharedref() type* buf)                                                        \
  {                                                                                  \
    GPUbarrier();                                                                    \
    ushort y = ll / N;                                                               \
    ushort x = ll % N;                                                               \
    Delta2 d = neighbors[x + offset];                                                \
    Delta dp = d.x;                                                                  \
    Delta dt = d.y;                                                                  \
    LOOP_UNROLL_ATTR for (unsigned int i = y; i < wgSize; i += (elems / N))          \
    {                                                                                \
      ChargePos readFrom = posBcast[i];                                              \
      uchar above = aboveThreshold[i];                                               \
      uint writeTo = N * i + x;                                                      \
      type v = null;                                                                 \
      if (expandFunc(above, x + offset)) {                                           \
        v = accessFunc(chargeMap, readFrom.gpad + dp, readFrom.time + dt);           \
      }                                                                              \
      buf[writeTo] = v;                                                              \
    }                                                                                \
    GPUbarrier();                                                                    \
  }
#else
#define DECL_FILL_SCRATCH_PAD(type, accessFunc)                                     \
  static GPUdi() void fillScratchPad_##type(                                        \
    GPUglobalref() const type* chargeMap,                                           \
    uint wgSize,                                                                    \
    uint,                                                                           \
    ushort ll,                                                                      \
    uint offset,                                                                    \
    uint N,                                                                         \
    GPUconstexprref() const Delta2* neighbors,                                      \
    GPUsharedref() const ChargePos* posBcast,                                       \
    GPUsharedref() type* buf)                                                       \
  {                                                                                 \
    if (ll >= wgSize) {                                                             \
      return;                                                                       \
    }                                                                               \
                                                                                    \
    ChargePos readFrom = posBcast[ll];                                              \
                                                                                    \
    GPUbarrier();                                                                   \
                                                                                    \
    for (unsigned int i = 0; i < N; i++) {                                          \
      Delta2 d = neighbors[i + offset];                                             \
      Delta dp = d.x;                                                               \
      Delta dt = d.y;                                                               \
                                                                                    \
      uint writeTo = N * ll + i;                                                    \
      buf[writeTo] = accessFunc(chargeMap, readFrom.gpad + dp, readFrom.time + dt); \
    }                                                                               \
                                                                                    \
    GPUbarrier();                                                                   \
  }

#define DECL_FILL_SCRATCH_PAD_COND(type, accessFunc, expandFunc, nameAppendix, null) \
  static GPUdi() void fillScratchPad##nameAppendix##_##type(                         \
    GPUglobalref() const type* chargeMap,                                            \
    uint wgSize,                                                                     \
    uint elems,                                                                      \
    ushort ll,                                                                       \
    uint offset,                                                                     \
    uint N,                                                                          \
    GPUconstexprref() const Delta2* neighbors,                                       \
    GPUsharedref() const ChargePos* posBcast,                                        \
    GPUsharedref() const uchar* aboveThreshold,                                      \
    GPUsharedref() type* buf)                                                        \
  {                                                                                  \
    if (ll >= wgSize) {                                                              \
      return;                                                                        \
    }                                                                                \
                                                                                     \
    ChargePos readFrom = posBcast[ll];                                               \
                                                                                     \
    GPUbarrier();                                                                    \
                                                                                     \
    for (unsigned int i = 0; i < N; i++) {                                           \
      Delta2 d = neighbors[i + offset];                                              \
      Delta dp = d.x;                                                                \
      Delta dt = d.y;                                                                \
                                                                                     \
      uint writeTo = N * ll + i;                                                     \
      uchar above = aboveThreshold[i];                                               \
      type v = null;                                                                 \
      if (expandFunc(above, i + offset)) {                                           \
        v = accessFunc(chargeMap, readFrom.gpad + dp, readFrom.time + dt);           \
      }                                                                              \
      buf[writeTo] = v;                                                              \
    }                                                                                \
                                                                                     \
    GPUbarrier();                                                                    \
  }
#endif

  DECL_FILL_SCRATCH_PAD(PackedCharge, CHARGE);
  DECL_FILL_SCRATCH_PAD(uchar, IS_PEAK);
  DECL_FILL_SCRATCH_PAD_COND(PackedCharge, CHARGE, innerAboveThreshold, Cond, PackedCharge(0));
  DECL_FILL_SCRATCH_PAD_COND(uchar, IS_PEAK, innerAboveThreshold, Cond, 0);
  DECL_FILL_SCRATCH_PAD_COND(PackedCharge, CHARGE, innerAboveThresholdInv, CondInv, PackedCharge(0));
  DECL_FILL_SCRATCH_PAD_COND(uchar, IS_PEAK, innerAboveThresholdInv, CondInv, 0);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
