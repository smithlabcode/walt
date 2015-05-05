#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include <queue>
#include <limits>

#include <fstream>
using std::ofstream;

#include "reference.hpp"

struct TEST_TIME {
  TEST_TIME() {
    get_region_start_t = 0;
    get_region_start_sum_time = 0;

    full_check_time_t = 0;
    full_check_sum_time = 0;

    num_of_full_check = 0;
  }
  clock_t get_region_start_t;
  uint64_t get_region_start_sum_time;

  clock_t full_check_time_t;
  uint64_t full_check_sum_time;

  uint64_t num_of_full_check;
};

struct BestMatch {
  BestMatch() {
    genome_pos = 0;
    times = 0;
    strand = '+';
    mismatch = MAX_INTEGER32;
  }
  BestMatch(const uint32_t& _genome_pos, const uint32_t& _times,
            const char& _strand, const uint32_t& _mismatch)
      : genome_pos(_genome_pos),
        times(_times),
        strand(_strand),
        mismatch(_mismatch) {
  }

  uint32_t genome_pos;
  uint32_t times;
  char strand;
  uint32_t mismatch;
};

struct CandidatePosition {
  CandidatePosition() {
    genome_pos = 0;
    strand = '+';
    mismatch = MAX_INTEGER32;
  }
  CandidatePosition(const uint32_t& _genome_pos, const char& _strand,
                    const uint32_t& _mismatch)
      : genome_pos(_genome_pos),
        strand(_strand),
        mismatch(_mismatch) {
  }

  bool operator<(const CandidatePosition & b) const {
    return mismatch < b.mismatch;
  }

  uint32_t genome_pos;
  char strand;
  uint32_t mismatch;
};

struct TopCandidates {
  TopCandidates() {
    size = 2;
  }

  TopCandidates(const uint32_t& _size)
      : size(_size) {
  }

  void SetSize(const uint32_t& _size) {
    size = _size;
  }

  bool Empty() {
    return candidates.empty();
  }

  void Clear() {
    while (!candidates.empty()) {
      candidates.pop();
    }
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

  CandidatePosition Top() {
    return candidates.top();
  }

  std::priority_queue<CandidatePosition> candidates;
  uint32_t size;
};

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      const HashTable& hash_table, BestMatch& best_match,
                      const uint32_t& seed_length, const char& strand,
                      TEST_TIME& test_time, const bool& AG_WILDCARD);

void PairEndMapping(const string& orginal_read, const Genome& genome,
                    const HashTable& hash_table, TopCandidates& top_match,
                    const uint32_t& max_mismatches, const uint32_t& seed_length,
                    const char& strand, const bool& AG_WILDCARD);

void MergePairedEndResults(
    const vector<vector<CandidatePosition> >& ranked_results,
    const vector<int>& ranked_results_size, const uint32_t& max_mismatches,
    const uint32_t& read_length, const uint32_t& frag_range,
    const Genome& genome, ofstream& fout);

#endif /* MAPPING_HPP_ */
