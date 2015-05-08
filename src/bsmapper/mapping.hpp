#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include <queue>
#include <limits>

#include <fstream>
using std::ofstream;

#include "reference.hpp"

/* BestMatch records the match position with the minimal number
 * of mismatches. times records how many genome positions with
 * this minimal number of mismatches. If times equals 1, means
 * the read is uniquely matched to the genome with the minimal
 * number of mismatches. */
struct BestMatch {
  BestMatch(const uint32_t& _genome_pos = 0, const uint32_t& _times = 0,
            const char& _strand = '+', const uint32_t& _mismatch = MAX_UINT32)
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
  TopCandidates(const uint32_t& _size = 100)
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

/* singled-end read mapping */
void SingleEndMapping(const string& org_read, const Genome& genome,
                      const HashTable& hash_table, const char& strand,
                      const bool& AG_WILDCARD, const uint32_t& seed_len,
                      BestMatch& best_match);

/* paired-end read mapping */
void PairEndMapping(const string& org_read, const Genome& genome,
                    const HashTable& hash_table, const char& strand,
                    const bool& AG_WILDCARD, const uint32_t& max_mismatches,
                    const uint32_t& seed_len, TopCandidates& top_match);

/* merge the mapping results from paired reads */
void MergePairedEndResults(
    const Genome& genome, const string& read_name, const string& read_seq1,
    const string& read_score1, const string& read_seq2,
    const string& read_score2,
    const vector<vector<CandidatePosition> >& ranked_results,
    const vector<int>& ranked_results_size, const uint32_t& read_len,
    const int& frag_range, const uint32_t& max_mismatches, ofstream& fout);

#endif /* MAPPING_HPP_ */
