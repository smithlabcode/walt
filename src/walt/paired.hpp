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
#include <vector>

/* CandidatePosition stores the candidate genome positions with number of
 * mismatches less or equal to max_mismatches */
struct CandidatePosition {
  CandidatePosition(const uint32_t gp = 0, const char st = '+',
                    const uint32_t mm = std::numeric_limits<uint32_t>::max()) :
    genome_pos(gp), strand(st), mismatch(mm) {}
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
  typedef std::priority_queue<CandidatePosition> cand_pq;

  TopCandidates(const uint32_t sz = 50) : max_size(sz) {}

  void SetSize(const uint32_t sz) {max_size = sz;}
  bool Empty() {return candidates.empty();}
  bool Full() {return candidates.size() >= max_size;}
  void Clear() {candidates = cand_pq();}
  void Pop() {candidates.pop();}
  CandidatePosition Top() {return candidates.top();}

  void Push(const CandidatePosition& cand) {
    if (candidates.size() < max_size)
      candidates.push(cand);
    else if (cand.mismatch < candidates.top().mismatch) {
      candidates.pop();
      candidates.push(cand);
    }
  }

  std::priority_queue<CandidatePosition> candidates;
  size_t max_size;
};

/* count the number of uniquely mapped, ambiguous mapped and
 * unmapped reads pairs */
struct StatPairedReads {
  StatPairedReads(const bool& _ambiguous, const bool& _unmapped,
                  const uint32_t& frag_range, const std::string &output_file,
                  const bool& SAM)
      : stat_single_reads_1(_ambiguous, _unmapped, output_file + "_1", SAM),
        stat_single_reads_2(_ambiguous, _unmapped, output_file + "_2", SAM) {
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
  std::vector<uint32_t> fragment_len_count;
};

/* paired-end read */
void ProcessPairedEndReads(const bool VERBOSE, const std::string &index_file,
                           const std::string &reads_file_p1,
                           const std::string &reads_file_p2,
                           const std::string &output_file,
                           const uint32_t& n_reads_to_process,
                           const uint32_t& max_mismatches, const uint32_t& b,
                           const std::string &adaptor,
                           const uint32_t& top_k, const int& frag_range,
                           const bool& ambiguous, const bool& unmapped,
                           const bool& SAM, const int& num_of_threads);

#endif /* PAIRED_HPP_ */
