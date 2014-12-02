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

/* for bisulfite sequence mapping, Cs are transfered to Ts*/
void C2T(const string& orginal_read, const uint32_t& read_len, string& read) {
  for (uint32_t i = 0; i < read_len; ++i) {
    if ('C' == orginal_read[i]) {
      read += 'T';
    } else {
      read += orginal_read[i];
    }
  }
}

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      BestMatch& best_match) {
  uint32_t read_len = orginal_read.size();
  if (read_len < HASHLEN)
    return;

  string read;
  C2T(orginal_read, read_len, read);

  uint32_t hash_value = getHashValue(&(read[0]));
  for (uint32_t i = 0; i < genome.size(); ++i) {
    const Chromosome& chrom = genome[i];
    HashTable::const_iterator it = chrom.hash_table.find(hash_value);
    if (it == chrom.hash_table.end())
      continue;
    for (uint32_t j = 0; j < it->second.size(); ++j) {
      if (it->second[j] + read_len >= chrom.length)
        return;

      /* check the position */
      uint32_t num_of_mismatch = 0;
      for (uint32_t q = it->second[j] + HASHLEN, p = HASHLEN; p < read_len;
          ++q, ++p) {
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
