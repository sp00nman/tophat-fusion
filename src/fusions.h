#ifndef FUSIONS_H
#define FUSIONS_H
/*
 *  fusions.h
 *  TopHat
 * 
 *  Adapted from junctions.h
 */

#include <cstdio>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <cstring>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>

#include "bwt_map.h"
using namespace std;

struct Fusion
{
  Fusion (uint32_t ref1, uint32_t ref2, uint32_t l, uint32_t r, uint32_t d = FUSION_FF)
  : refid1(ref1), refid2(ref2), left(l), right(r), dir(d)
  {}
  
Fusion() : refid1(0), refid2(0), left(0), right(0), dir(FUSION_FF)
  {}
  
  uint32_t refid1;
  uint32_t refid2;
  uint32_t left;
  uint32_t right;
  uint32_t dir;
  
  bool operator<(const Fusion& rhs) const
  {
    if (refid1 < rhs.refid1)
      return true;
    else if (refid1 > rhs.refid1)
      return false;
    
    if (refid2 < rhs.refid2)
      return true;
    else if (refid2 > rhs.refid2)
      return false;
    
    if (left < rhs.left)
      return true;
    else if (left > rhs.left)
      return false;

    if (right < rhs.right)
      return true;
    else if (right > rhs.right)
      return false;

    if (dir < rhs.dir)
      return true;
    else if(dir > rhs.dir)
      return false;

    return false;
  }
	
  bool operator==(const Fusion& rhs) const
  {
    return  refid1 == rhs.refid1 &&
    refid2 == rhs.refid2 &&
    left == rhs.left &&
    right == rhs.right &&
    dir == rhs.dir;
  }
};

struct FusionPairSupport
{
FusionPairSupport(int l, int r) :
  ldist(l),
    rdist(r)
  {}

  bool operator<(const FusionPairSupport& rhs) const
  {
    if (abs(ldist) + abs(rdist) < abs(rhs.ldist) + abs(rhs.rdist))
      return true;
    else
      return false;
  }

  int ldist;
  int rdist;
};

struct FusionSimpleStat
{
  uint32_t count; // # of reads that support the fusion
  uint32_t edit_dist; // the smallest edit dist among the supporting reads
  bool skip; //
};

struct FusionStat
{
  uint32_t count; // # of reads that support the fusion
  uint32_t unsupport_count; // # of reads that unsupport the fusion, that is, reads that span one chromosome.
  uint32_t unsupport_count_pair;
  uint32_t pair_count; // # of pairs that support the fusion
  uint32_t pair_count_fusion; // # of pairs that support the fusion where one read spans the fusion
  float symm;

  string chr1_seq;
  string chr2_seq;
  
  static const uint32_t NUM_BASES = 50;
  uint32_t left_bases[NUM_BASES];
  uint32_t right_bases[NUM_BASES];

  uint32_t diffs[5];

  uint32_t left_ext;
  uint32_t right_ext;

  vector<FusionPairSupport> vPairSupport;

  /*
   * daehwan - this is a metadata indicating whether Fusion is reversed ..
   * it's true that this is not a good way to ...
   */
  bool reversed;

  FusionStat()
  {
    count = 0;
    unsupport_count = 0;
    unsupport_count_pair = 0;
    pair_count = 0;
    pair_count_fusion = 0;
    symm = 0.0f;
    memset(left_bases, 0, sizeof(left_bases));
    memset(right_bases, 0, sizeof(right_bases));
    memset(diffs, 0, sizeof(diffs));

    left_ext = right_ext = 0;

    reversed = false;
  }
};

struct fusion_comparison 
{
  bool operator()(const Fusion& lhs, const Fusion& rhs)
  {
    if (lhs.left != rhs.left)
      return lhs.left < rhs.left;

    if (lhs.right != rhs.right)
      return lhs.right < rhs.right;

    if (lhs.refid1 != rhs.refid1)
      return lhs.refid1 < rhs.refid1;

    if (lhs.refid2 != rhs.refid2)
      return lhs.refid2 < rhs.refid2;

    if (lhs.dir != rhs.dir)
      return lhs.dir < rhs.dir;

    return false;
  }
};

// this is used in segment_juncs
typedef std::map<Fusion, FusionSimpleStat> FusionSimpleSet;

// this is used in tophat_reports
typedef std::map<Fusion, FusionStat, fusion_comparison> FusionSet;

void fusions_from_alignment(const BowtieHit& bh, FusionSet& fusions, RefSequenceTable& rt, bool update_stat = false);
void unsupport_fusions(const BowtieHit& bh, FusionSet& fusions, FusionSet& fusions_ref);
void print_fusions(FILE* fusions_out, FusionSet& fusions, RefSequenceTable& ref_sequences);
void fusions_from_spliced_hit(const BowtieHit& bh, vector<Fusion>& fusions, bool auto_sort = true);
void pair_support(const HitsForRead& left_hits, const HitsForRead& right_hits, FusionSet& fusions, FusionSet& fusions_ref);


#endif
