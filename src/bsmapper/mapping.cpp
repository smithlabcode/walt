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
    if ('N' == orginal_read[i]) {
      read += getNT(3);  // in rampbs N set to 3.
    } else if ('C' == orginal_read[i]) {
      read += 'T';
    } else {
      read += orginal_read[i];
    }
  }
}

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      const HashTable& hash_table, BestMatch& best_match,
                      const uint32_t& HASHLEN) {
  uint32_t read_len = orginal_read.size();
  if (read_len < HASHLEN)
    return;

  string read;
  C2T(orginal_read, read_len, read);

  uint64_t hash_value = getHashValue(&(read[0]), HASHLEN);
  HashTable::const_iterator it = hash_table.find(hash_value);
  if (it == hash_table.end())
    return;

  for (uint32_t j = 0; j < it->second.size(); ++j) {
    const Chromosome& chrom = genome[it->second[j].first];
    if (it->second[j].second + read_len >= chrom.length)
      return;

    /* check the position */
    uint32_t num_of_mismatch = 0;
    for (uint32_t q = it->second[j].second + HASHLEN, p = HASHLEN; p < read_len;
        ++q, ++p) {
      if (chrom.sequence[q] != read[p]) {
        num_of_mismatch++;
      }
      if (num_of_mismatch > best_match.mismatch)
        break;
    }

    if (num_of_mismatch < best_match.mismatch) {
      best_match = BestMatch(it->second[j].first, it->second[j].second, 1,
                             num_of_mismatch);
    } else if (best_match.mismatch == num_of_mismatch) {
      best_match.times++;
    }
  }
}
