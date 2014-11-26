#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include <limits>
#include "reference.hpp"

struct BestMatch {
  BestMatch() {
    times = 0;
    mismatch = std::numeric_limits<uint32_t>::max();
  }
  BestMatch(const GenomePosition& _genome_pos, const uint32_t& _times,
            const uint32_t& _mismatch)
      : genome_pos(_genome_pos),
        times(_times),
        mismatch(_mismatch) {
  }
//  BestMatch& operator=(const BestMatch& other) {
//    this->diag = other.diag;
//    this->times = other.times;
//    this->mismatch = other.mismatch;
//    this->start_pos = other.start_pos;
//    this->strand = other.strand;
//    return *this;
//  }

  bool operator<(const BestMatch& other) {
    return this->mismatch < other.mismatch;
  }

  GenomePosition genome_pos;
  uint32_t times;
  uint32_t mismatch;
};

struct DiagSize {
  DiagSize(const uint32_t& _diag, const uint32_t& _chrom_id,
           const uint32_t& _size)
      : diag(_diag),
        chrom_id(_chrom_id),
        size(_size) {
  }
  uint32_t diag;
  uint32_t chrom_id;
  uint32_t size;
  bool operator<(const DiagSize& b) const {
    return size > b.size;
  }
};

void SingleEndMapping(const char* read, const Genome& genome,
                      const int& num_top_diags);

#endif /* MAPPING_HPP_ */
