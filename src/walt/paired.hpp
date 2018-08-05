/*
 *    The head file for mapping paired-end reads
 *
 *    Copyright (C) 2015 University of Southern California
 *                       Andrew D. Smith and Ting Chen
 *
 *    Authors: Haifeng Chen, Andrew D. Smith and Ting Chen
 *
 *    This file is part of WALT.
 *
 *    WALT is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    WALT is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with WALT.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PAIRED_HPP_
#define PAIRED_HPP_

#include "mapping.hpp"

#include <queue>

/* CandidatePosition stores the candidate genome positions with number of
 * mismatches less or equal to max_mismatches */
struct CandidatePosition {
  CandidatePosition(const uint32_t& _genome_pos = 0, const char& _strand = '+',
                    const uint32_t& _mismatch = MAX_UINT32)
      : genome_pos(_genome_pos),
        strand(_strand),
        mismatch(_mismatch) {
  }

  bool operator<(const CandidatePosition& b) const {
    return mismatch < b.mismatch;
  }

  uint32_t genome_pos;
  char strand;
  uint32_t mismatch;
};

/* TopCandidates is a priority_queue which stores the top-k candidate
 * positions. When mapping paired-end reads, for each read in the pair,
 * the top-k positions (with minimal mismatches) are recorded. Then using
 * the top-k positions in each of them to find the best pair match. */
struct TopCandidates {
  TopCandidates(const uint32_t& _size = 50)
      : size(_size) {
  }

  void SetSize(const uint32_t& _size) {
    size = _size;
  }

  bool Empty() {
    return candidates.empty();
  }

  bool Full() {
    return candidates.size() >= size;
  }

  void Clear() {
    while (!candidates.empty()) {
      candidates.pop();
    }
  }

  CandidatePosition Top() {
    return candidates.top();
  }

  void Push(const CandidatePosition& cand) {
    if (candidates.size() < size) {
      candidates.push(cand);
    } else {
      if (cand.mismatch < candidates.top().mismatch) {
        candidates.pop();
        candidates.push(cand);
      }
    }
  }

  void Pop() {
    candidates.pop();
  }

  std::priority_queue<CandidatePosition> candidates;
  uint32_t size;
};

/* count the number of uniquely mapped, ambiguous mapped and
 * unmapped reads pairs */
struct StatPairedReads {
  StatPairedReads(const bool& _ambiguous, const bool& _unmapped,
                  const uint32_t& frag_range,
                  const bool& SAM)
      : stat_single_reads_1(_ambiguous, _unmapped, SAM),
        stat_single_reads_2(_ambiguous, _unmapped, SAM) {
    total_read_pairs = 0;
    unique_mapped_pairs = 0;
    ambiguous_mapped_pairs = 0;
    unmapped_pairs = 0;

    fragment_len_count.resize(frag_range + 1);
    for (uint32_t i = 0; i < fragment_len_count.size(); ++i) {
      fragment_len_count[i] = 0;
    }
  }

  uint32_t total_read_pairs;
  uint32_t unique_mapped_pairs;
  uint32_t ambiguous_mapped_pairs;
  uint32_t unmapped_pairs;

  StatSingleReads stat_single_reads_1;
  StatSingleReads stat_single_reads_2;

  // distribution of fragment length
  vector<uint32_t> fragment_len_count;
};

/* paired-end read */
void ProcessPairedEndReads(const string& index_file,
                           const string& reads_file_p1,
                           const string& reads_file_p2,
                           FILE * fout, FILE * fstat,
                           const uint32_t& n_reads_to_process,
                           const uint32_t& max_mismatches, const uint32_t& b,
                           const string& adaptor, const bool& PBAT,
                           const uint32_t& top_k, const int& frag_range,
                           FILE * fambiguous, FILE * funmapped,
                           const bool& SAM);

#endif /* PAIRED_HPP_ */
