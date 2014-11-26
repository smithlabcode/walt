#include "mapping.hpp"

#include <queue>
using std::priority_queue;

string ReverseComplimentStrand(const string& read) {
  string reverse_complement_read;
  uint32_t read_len = read.size();
  for (uint32_t i = 0; i < read_len; ++i) {
    reverse_complement_read += complimentBase(read[read_len - i - 1]);
  }
  return reverse_complement_read;
}

void SingleEndMapping(const char* read, const Genome& genome,
                      BestMatch& best_match) {
  uint32_t read_len = strlen(read);
  if (read_len < HASHLEN)
    return;

  uint32_t hash_value = getHashValue(&(read[0]));
  for (uint32_t i = 0; i < genome.size(); ++i) {
    const Chromosome& chrom = genome[i];
    HashTable::const_iterator it = chrom.hash_table.find(hash_value);
    for (uint32_t j = 0; j < it->second.size(); ++j) {
      if (it->second[j] + read_len > chrom.length)
        return;

      /* check the position */
      uint32_t num_of_mismatch = 0;
      for (uint32_t q = it->second[j], p = 0; p < read_len; ++q, ++p) {
        if (chrom.sequence[q] != read[p]) {
          num_of_mismatch++;
        }
        if (num_of_mismatch > best_match.mismatch)
          break;
      }

      if (num_of_mismatch < best_match.mismatch) {
        best_match = BestMatch(i, it->second[j], 1, num_of_mismatch);
      } else if (best_match.mismatch == num_of_mismatch) {
        best_match.times++;
      }
    }
  }
}
