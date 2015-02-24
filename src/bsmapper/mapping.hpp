#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include <limits>
#include "reference.hpp"

struct BestMatch {
  BestMatch() {
    chrom_id = 0;
    chrom_pos = 0;
    times = 0;
    mismatch = std::numeric_limits<uint32_t>::max();
  }
  BestMatch(const uint32_t&_chrom_id, const uint32_t& _chrom_pos,
            const uint32_t& _times, const uint32_t& _mismatch)
      : chrom_id(_chrom_id),
        chrom_pos(_chrom_pos),
        times(_times),
        mismatch(_mismatch) {
  }

  uint32_t chrom_id;
  uint32_t chrom_pos;
  uint32_t times;
  uint32_t mismatch;
};

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      const HashTable& hash_table, BestMatch& best_match,
                      const uint32_t& seed_length);

#endif /* MAPPING_HPP_ */
