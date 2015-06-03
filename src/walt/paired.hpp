/*
 * paired-end read mapping
 * */
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
    return candidates.size() == size;
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
                  const string& output_file, const bool& SAM)
      : stat_single_reads_1(_ambiguous, _unmapped, output_file + "_1", SAM),
        stat_single_reads_2(_ambiguous, _unmapped, output_file + "_2", SAM) {
    total_read_pairs = 0;
    unique_mapped_pairs = 0;
    ambiguous_mapped_pairs = 0;
    unmapped_pairs = 0;
  }

  uint32_t total_read_pairs;
  uint32_t unique_mapped_pairs;
  uint32_t ambiguous_mapped_pairs;
  uint32_t unmapped_pairs;

  StatSingleReads stat_single_reads_1;
  StatSingleReads stat_single_reads_2;
};

/* paired-end read */
void ProcessPairedEndReads(const string& command, const string& index_file,
                           const string& reads_file_p1,
                           const string& reads_file_p2,
                           const string& output_file,
                           const uint32_t& n_reads_to_process,
                           const uint32_t& max_mismatches,
                           const string& adaptor, const uint32_t& top_k,
                           const int& frag_range, const bool& ambiguous,
                           const bool& unmapped, const bool& SAM);

#endif /* PAIRED_HPP_ */
