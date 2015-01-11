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

void Validation(
    const vector<string>& c2t_read_seqs, const Genome& genome,
    const HashTable& hash_table,
    const unordered_map<uint64_t, vector<uint32_t> >& hash_value_map2_reads,
    const uint32_t& seed_i, vector<BestMatch>& map_results) {
  for (unordered_map<uint64_t, vector<uint32_t> >::const_iterator it =
      hash_value_map2_reads.begin(); it != hash_value_map2_reads.end(); ++it) {
    HashTable::const_iterator it2 = hash_table.find(it->first);
    for (uint32_t j = 0; j < it2->second.size(); ++j) {
      for (uint32_t i = 0; i < it->second.size(); ++i) {
        const GenomePosition& gp = it2->second[j];
        const Chromosome& chrom = genome[gp.chrom_id];
        const uint32_t& read_id = it->second[i];
        if (gp.chrom_pos >= seed_i) {
          uint32_t chrom_pos = gp.chrom_pos - seed_i;
          uint32_t read_len = c2t_read_seqs[read_id].size();
          if (chrom_pos + read_len < chrom.length) {
            /* check the position */
            uint32_t num_of_mismatch = 0;
            for (uint32_t q = chrom_pos, p = 0;
                p < read_len && num_of_mismatch <= map_results[read_id].mismatch;
                ++q, ++p) {
              if (chrom.sequence[q] != c2t_read_seqs[it->second[i]][p]) {
                num_of_mismatch++;
              }
            }
            if (num_of_mismatch < map_results[read_id].mismatch) {
              map_results[read_id] = BestMatch(gp.chrom_id, chrom_pos, 1,
                                               num_of_mismatch);
            } else if (map_results[read_id].mismatch == num_of_mismatch) {
              if (map_results[read_id].chrom_id != gp.chrom_id
                  || map_results[read_id].chrom_pos != chrom_pos) {
                map_results[read_id].chrom_id = gp.chrom_id;
                map_results[read_id].chrom_pos = chrom_pos;
                map_results[read_id].times++;
              }
            }
          }
        }
      }
    }
  }
}

void MappingAllReads(const vector<string>& read_seqs,
                     const uint32_t& num_of_reads, const Genome& genome,
                     const HashTable& hash_table,
                     vector<BestMatch>& map_results) {
  vector<string> c2t_read_seqs(num_of_reads);
  for (uint32_t i = 0; i < num_of_reads; ++i) {
    C2T(read_seqs[i], read_seqs[i].size(), c2t_read_seqs[i]);
    //cerr << c2t_read_seqs[i] << endl;
  }

  for (uint32_t seed_i = 0; seed_i < 7; ++seed_i) {
    //cerr << seed_i << " 7 " << endl;
    unordered_map<uint64_t, vector<uint32_t> > hash_value_map2_reads;
    for (uint32_t i = 0; i < num_of_reads; ++i) {
      string read_seed = c2t_read_seqs[i].substr(seed_i);
      uint64_t hash_value = getHashValue(read_seed.c_str());
      HashTable::const_iterator it = hash_table.find(hash_value);
      if (it == hash_table.end())
        continue;
      hash_value_map2_reads[hash_value].push_back(i);
    }
    Validation(c2t_read_seqs, genome, hash_table, hash_value_map2_reads, seed_i,
               map_results);
  }
}
