#include <utility>
#include <string>
#include <vector>

#include "mapping.hpp"

#include <tr1/unordered_set>
using std::tr1::unordered_set;

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

int smaller(const GenomePosition& a, const uint32_t& seed_i,
            const GenomePosition& b, const uint32_t& seed_j) {
  if (a.chrom_id < b.chrom_id)
    return 1;
  else if (a.chrom_id > b.chrom_id)
    return -1;
  else {
    if (a.chrom_pos - seed_i < b.chrom_pos - seed_j - HASHLEN)
      return 1;
    else if (a.chrom_pos - seed_i > b.chrom_pos - seed_j - HASHLEN)
      return -1;
    return 0;
  }
}

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      const HashTable& hash_table, BestMatch& best_match) {
  uint32_t read_len = orginal_read.size();
  if (read_len < HASHLEN2)
    return;

  string read;
  C2T(orginal_read, read_len, read);
  vector<pair<uint32_t, uint32_t> > hash_values_1;
  for (uint32_t seed_i = 0; seed_i < 7; ++seed_i) {
    string read_seed = read.substr(seed_i);
    uint32_t hash_value = getHashValue(read_seed.c_str());
    HashTable::const_iterator it = hash_table.find(hash_value);
    if (it == hash_table.end())
      continue;
    hash_values_1.push_back(make_pair(hash_value, seed_i));
  }

  vector<pair<uint32_t, uint32_t> > hash_values_2;
  for (uint32_t seed_i = 0; seed_i < 7; ++seed_i) {
    string read_seed = read.substr(seed_i + HASHLEN);
    uint32_t hash_value = getHashValue(read_seed.c_str());
    HashTable::const_iterator it2 = hash_table.find(hash_value);
    if (it2 == hash_table.end())
      continue;
    hash_values_2.push_back(make_pair(hash_value, seed_i));
  }

  for (uint32_t i = 0; i < hash_values_1.size(); ++i) {
    for (uint32_t j = 0; j < hash_values_2.size(); ++j) {
      HashTable::const_iterator it1 = hash_table.find(hash_values_1[i].first);
      HashTable::const_iterator it2 = hash_table.find(hash_values_2[j].first);
      uint32_t p = 0, q = 0;
      while (p < it1->second.size() && q < it2->second.size()) {
        int s = smaller(it1->second[p], hash_values_1[i].second, it2->second[q],
                        hash_values_2[j].second);
        if (s == 1) {
          p++;
        } else if (s == -1) {
          q++;
        } else {
          const GenomePosition& gp = it1->second[p];
          const Chromosome& chrom = genome[gp.chrom_id];
          uint32_t chrom_pos = gp.chrom_pos - hash_values_1[i].second;
          if (chrom_pos + read_len < chrom.length) {
            /* check the position */
            uint32_t num_of_mismatch = 0;
            for (uint32_t q = chrom_pos, p = 0;
                p < read_len && num_of_mismatch <= best_match.mismatch;
                ++q, ++p) {
              if (chrom.sequence[q] != read[p]) {
                num_of_mismatch++;
              }
            }

            if (num_of_mismatch < best_match.mismatch) {
              best_match = BestMatch(gp.chrom_id, chrom_pos, 1,
                                     num_of_mismatch);
            } else if (best_match.mismatch == num_of_mismatch) {
              if (best_match.chrom_id != gp.chrom_id
                  || best_match.chrom_pos != chrom_pos) {
                best_match.chrom_id = gp.chrom_id;
                best_match.chrom_pos = chrom_pos;
                best_match.times++;
              }
            }
          }
        }
        p++;
        q++;
      }
    }
  }
}

