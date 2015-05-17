/*
 * paired-end read mapping
 * */
#ifndef PAIRED_HPP_
#define PAIRED_HPP_

#include <algorithm>

#include "mapping.hpp"

/* CandidatePosition stores the candidate genome positions with number of
 * mismatches less or equal to max_mismatches */
struct CandidatePosition {
  CandidatePosition(const uint32_t& _genome_pos = 0, const char& _strand = '+',
                    const uint32_t& _mismatch = MAX_UINT32)
      : genome_pos(_genome_pos),
        strand(_strand),
        mismatch(_mismatch) {
  }
  void Output() const {
    fprintf(stderr, "%u %c %u\n", genome_pos, strand, mismatch);
  }

  uint32_t genome_pos;
  char strand;
  uint32_t mismatch;
};

/* TopCandidates is a priority_queue which stores the top-k candidate
 * positions. When mapping paired-end reads, for each read in the pair,
 * the top-k positions (with minimal mismatches) are recorded. Then using
 * the top-k positions in each of them to find the best pair match. */

struct SORTCMP {
  bool operator()(const CandidatePosition& a, const CandidatePosition& b) {
    return a.mismatch < b.mismatch;
  }
};

struct TopCandidates_Heap {
  TopCandidates_Heap(const uint32_t& _size = 100) {
    candidates.resize(_size);
    capacity = _size;
  }

  void Clear() {
    cur_size = 0;
  }

  CandidatePosition Top() {
    std::make_heap(candidates.begin(), candidates.begin() + cur_size,
                   SORTCMP());

    return candidates.front();
  }

  void Push(const CandidatePosition& cand) {
    if (cur_size < capacity) {
      candidates[cur_size] = cand;
      cur_size++;
    } else {
      if (cand.mismatch < Top().mismatch) {
        std::pop_heap(candidates.begin(), candidates.begin() + cur_size,
                      SORTCMP());
        candidates[capacity - 1] = cand;
      }
    }
  }

  void Sort() {
    std::make_heap(candidates.begin(), candidates.begin() + cur_size,
                   SORTCMP());
    std::sort_heap(candidates.begin(), candidates.begin() + cur_size,
                   SORTCMP());
  }

  void Output() {
    for (uint32_t i = 0; i < cur_size; ++i) {
      candidates[i].Output();
    }
  }

  vector<CandidatePosition> candidates;
  uint32_t cur_size;
  uint32_t capacity;
};

/* paired-end read */
void ProcessPairedEndReads(const string& index_file,
                           const string& reads_file_p1,
                           const string& reads_file_p2,
                           const string& output_file,
                           const uint32_t& n_reads_to_process,
                           const uint32_t& max_mismatches,
                           const uint32_t& read_len, const uint32_t& seed_len,
                           const uint32_t& top_k, const int& frag_range);

#endif /* PAIRED_HPP_ */
