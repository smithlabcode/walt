#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include <limits>
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
    mismatch = std::numeric_limits < uint32_t > ::max();
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

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      const HashTable& hash_table, BestMatch& best_match,
                      const uint32_t& seed_length, const char& strand,
                      TEST_TIME& test_time);

#endif /* MAPPING_HPP_ */
