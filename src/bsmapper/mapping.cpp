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

uint32_t GetDiag(const uint32_t& chrom_pos, const uint32_t& seed_pos,
                 const uint32_t& read_len) {
  uint32_t diag = chrom_pos + read_len - seed_pos;

  return diag;
}

void GetDiagsSize(const vector<uint32_t>& hash_values, const Chromosome& chrom,
                  const uint32_t& read_len,
                  unordered_map<uint32_t, uint32_t>& diags_size) {
  /* count how many matched seeds in each diag */
  for (uint32_t i = 0; i < hash_values.size(); ++i) {
    HashTable::const_iterator it = chrom.hash_table.find(hash_values[i]);
    if (it == chrom.hash_table.end())
      continue;

    /* ignore large bucket*/
    if (it->second.size() > 5000)
      continue;

    for (uint32_t val = 0; val < it->second.size(); ++val) {
      diags_size[GetDiag(it->second[val], i, read_len)]++;
    }
  }
}

void GetTopDiags(const unordered_map<uint32_t, uint32_t>& diags_size,
                 priority_queue<DiagSize>& top_diags_queue,
                 const uint32_t& chrom_id, const int& num_top_diags) {
  /* select the top diags */
  uint32_t num_top_diags_threshold = num_top_diags + 1;

  for (unordered_map<uint32_t, uint32_t>::const_iterator it =
      diags_size.begin(); it != diags_size.end(); it++) {
    top_diags_queue.push(DiagSize(it->first, chrom_id, it->second));
    if (top_diags_queue.size() == num_top_diags_threshold) {
      top_diags_queue.pop();
    }
  }
}

void SingleEndMapping(const char* read, const Genome& genome,
                      const int& num_top_diags) {
  uint32_t read_len = strlen(read);
  uint32_t num_of_seeds = read_len - HASHLEN + 1;
  vector<uint32_t> hash_values(num_of_seeds);
  for (uint32_t i = 0; i < num_of_seeds; ++i) {
    hash_values[i] = getHashValue(&(read[i]));
  }

  priority_queue<DiagSize> top_diags_queue;
  for (uint32_t i = 0; i < genome.size(); ++i) {
    unordered_map<uint32_t, uint32_t> diags_size;
    GetDiagsSize(hash_values, genome[i], read_len, diags_size);
    GetTopDiags(diags_size, top_diags_queue, i, num_top_diags);
  }

  vector<DiagSize> top_diags;
  while (!top_diags_queue.empty()) {
    top_diags.push_back(top_diags_queue.top());
    top_diags_queue.pop();
  }

  /* select the best match */
  BestMatch best_match;

  int top_diags_size = top_diags.size();
  for (int i = top_diags_size - 1; i >= 0; --i) {
    uint32_t chrom_id = top_diags[i].chrom_id;
    uint32_t diag = top_diags[i].diag;
    uint32_t num_of_seed_match = top_diags[i].size;

    if (diag >= genome[chrom_id].length)
      continue;

    if (diag < read_len)
      continue;

    uint32_t chrom_pos = diag - read_len;

    GenomePosition genome_pos = make_pair(chrom_id, chrom_pos);

    if (num_of_seed_match == num_of_seeds) {
      /* exact match */
      if (best_match.mismatch > 0) {
        best_match = BestMatch(genome_pos, 1, 0);
        continue;
      }

      if (best_match.mismatch == 0) {
        cerr << "[exact match to multiple positions, ignore this read]" << endl;
        return;
      }
    }

    if (best_match.mismatch == 0)
      continue;

    /* mismatch */

    uint32_t num_of_mismatch = 0;
    uint32_t k = chrom_pos;
    for (uint32_t j = 0; j < read_len; ++j) {
      if (genome[chrom_id].sequence[k] != read[j]) {
        num_of_mismatch++;
      }
      if (num_of_mismatch > best_match.mismatch)
        break;
      k++;
    }

    if (num_of_mismatch < best_match.mismatch) {
      best_match = BestMatch(genome_pos, 1, num_of_mismatch);
    } else if (best_match.mismatch == num_of_mismatch) {
      best_match.times++;
    }
  }
  if (best_match.times > 1) {
    cerr << "[best match to multiple positions, ignore this read]" << endl;
    return;
  }
}
