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

GenomePosition GetDiag(const uint32_t& chrom_pos, const uint32_t& seed_pos,
                       const uint32_t& chrom_id, const uint32_t& read_len) {
  GenomePosition pos;
  uint32_t diag = chrom_pos + read_len - seed_pos;

  return make_pair(chrom_id, diag);
}

uint32_t GetGenomeStartPos(const uint32_t& diag, const uint32_t& read_len) {
  if (read_len > diag)
    return 0;

  return diag - read_len;
}

void GetDiagsSize(const char* read, uint32_t& read_len,
                  unordered_map<GenomePosition, uint32_t>& diags_size,
                  const Genome& genome) {
  /* count how many matched seeds in each diag */
  uint32_t hash_value = 0, num_of_seeds = read_len - HASHLEN + 1;
  for (uint32_t i = 0; i < num_of_seeds; ++i) {
    hash_value = getHashValue(&(read[i]));
    for (uint32_t j = 0; j < genome.size(); ++j) {
      HashTable::const_iterator it = genome[j].hash_table.find(hash_value);
      if (it == genome[j].hash_table.end())
        continue;

      /* ignore large bucket*/
      if (it->second.size() > 5000)
        continue;

      for (uint32_t val = 0; val < it->second.size(); ++val) {
        diags_size[GetDiag(it->second[val], i, j, read_len)]++;
      }
    }
  }
}

void GetTopDiags(unordered_map<GenomePosition, uint32_t>& diags_size,
                 vector<GenomePosition>& top_diags, const int& num_top_diags) {
  /* select the top diags */
  priority_queue<DiagSize> top_diags_queue;
  uint32_t num_top_diags_threshold = num_top_diags + 1;

  for (unordered_map<GenomePosition, uint32_t>::iterator it =
      diags_size.begin(); it != diags_size.end(); it++) {
    top_diags_queue.push(DiagSize(it->first, it->second));
    if (top_diags_queue.size() == num_top_diags_threshold) {
      top_diags_queue.pop();
    }
  }

  while (!top_diags_queue.empty()) {
    top_diags.push_back(top_diags_queue.top().diag);
    top_diags_queue.pop();
  }
}

void SingleEndMapping(const string& read, const Genome& genome,
                      const int& num_top_diags) {
  uint32_t read_len = read.size();
  uint32_t num_of_seeds = read_len - HASHLEN + 1;

  unordered_map<GenomePosition, uint32_t> diags_size;
  GetDiagsSize(read.c_str(), read_len, diags_size, genome);

  vector<GenomePosition> top_diags;
  GetTopDiags(diags_size, top_diags, num_top_diags);

  /* select the best match */
  BestMatch best_match;

  int top_diags_size = top_diags.size();
  for (int i = top_diags_size - 1; i >= 0; --i) {
    uint32_t chrom_id = top_diags[i].first;
    uint32_t chrom_pos = top_diags[i].second;

    if (chrom_pos >= genome[chrom_id].length)
      continue;

    uint32_t start_pos = GetGenomeStartPos(chrom_pos, read_len);

    GenomePosition genome_pos = make_pair(chrom_id, start_pos);

    if (diags_size[top_diags[i]] == num_of_seeds) {
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
    uint32_t k = start_pos;
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
