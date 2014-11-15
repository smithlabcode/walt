#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include "util.hpp"
#include "index.hpp"
#include "reference.hpp"

struct DiagSize {
  DiagSize(const uint32_t& _diag, const uint32_t& _size)
      : diag(_diag),
        size(_size) {
  }
  uint32_t diag;
  uint32_t size;
  bool operator<(const DiagSize& b) const {
    return size > b.size;
  }
};


/* load reads from reads file, each time load n_reads_to_process reads,
 * start from  read_start_idx */
void LoadReadsFromFastqFile(const string &filename,
                            const uint64_t read_start_idx,
                            const uint64_t n_reads_to_process,
                            vector<string>& read_names,
                            vector<string>& read_seqs);

/* map the read to the genome */
void mapping(const Genome* genome, const HashTable* hash_table,
             const char* read, const int& num_top_diags);

#endif /* MAPPING_HPP_ */
