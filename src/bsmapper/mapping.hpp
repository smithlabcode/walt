#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include "util.hpp"
#include "index.hpp"
#include "reference.hpp"

#include <tr1/unordered_map>
#include <limits>

using std::make_pair;
using std::tr1::unordered_map;

struct BestMatch {
  BestMatch() {
    diag = 0;
    times = 0;
    mismatch = std::numeric_limits<uint32_t >::max();
    start_pos = std::numeric_limits<uint32_t>::max();
    strand = '+';
  }
  BestMatch(const uint32_t& _diag, const uint32_t& _times,
            const uint32_t& _mismatch, const uint32_t& _start_pos,
            char& _strand)
      : diag(_diag),
        times(_times),
        mismatch(_mismatch),
        start_pos(_start_pos),
        strand(_strand) {
  }
  BestMatch& operator=(const BestMatch& other) {
    this->diag = other.diag;
    this->times = other.times;
    this->mismatch = other.mismatch;
    this->start_pos = other.start_pos;
    this->strand = other.strand;
    return *this;
  }

  bool operator<(const BestMatch& other) {
    return this->mismatch < other.mismatch;
  }

  uint32_t diag;
  uint32_t times;
  uint32_t mismatch;
  uint32_t start_pos;
  char strand;
};

struct DiagSize {
  DiagSize(const uint32_t& _diag, const uint32_t& _size, const char& _strand)
      : diag(_diag),
        size(_size),
        strand(_strand) {
  }
  uint32_t diag;
  uint32_t size;
  char strand;
  bool operator<(const DiagSize& b) const {
    return size > b.size;
  }
};

class Mapping {
 public:
  Mapping(const Genome* _genome, const HashTable* _hash_table,
          const int& _num_top_diags)
      : genome(_genome),
        hash_table(_hash_table),
        num_top_diags(_num_top_diags) {
  }
  /* map the read to the genome */
  void SingleEndMapping(const string& read);

 private:
  void SingleEndMapping(const char* read);
  string ReverseComplimentStrand(const string& read);

  uint32_t GetDiag(const uint32_t& genome_pos, const uint32_t& seed_pos,
                   const uint32_t& read_len);
  uint32_t GetGenomeStartPos(const uint32_t& diag, const uint32_t& read_len);
  void GetDiagsSize(const char* read, uint32_t& read_len,
                    unordered_map<uint32_t, uint32_t>& diags_size);
  void GetTopDiags(const unordered_map<uint32_t, uint32_t>& diags_size_pos,
                   const unordered_map<uint32_t, uint32_t>& diags_size_neg,
                   vector<pair<uint32_t, char> >& top_diags);

  const Genome* genome;
  const HashTable* hash_table;
  int num_top_diags;

};

#endif /* MAPPING_HPP_ */
