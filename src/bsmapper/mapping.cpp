#include "mapping.hpp"

#include <queue>
#include <tr1/unordered_map>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::tr1::unordered_map;
using std::priority_queue;

uint32_t Mapping::GetDiag(const uint32_t& genome_pos, const uint32_t& seed_pos,
                          const uint32_t& read_len) {
  return genome_pos + read_len - seed_pos;
}

uint32_t Mapping::GetGenomeStartPos(const uint32_t& diag,
                                    const uint32_t& read_len) {
  if (read_len > diag)
    return 0;
  return diag - read_len;
}

void Mapping::GetDiagsSize(const char* read, uint32_t& read_len,
                           unordered_map<uint32_t, uint32_t>& diags_size) {
  /* count how many matched seeds in each diag */
  uint32_t hash_value = 0, num_of_seeds = read_len - HASHLEN + 1;
  for (uint32_t i = 0; i < num_of_seeds; ++i) {
    hash_value = getHashValue(&(read[i]));
    HashTable::const_iterator it = hash_table->find(hash_value);
    if (it == hash_table->end())
      continue;

    /* ignore large bucket*/
    if (it->second.size() > 5000)
      continue;

    for (uint32_t val = 0; val < it->second.size(); ++val) {
      diags_size[GetDiag(it->second[val], i, read_len)]++;
    }
  }
}

void Mapping::GetTopDiags(
    const unordered_map<uint32_t, uint32_t>& diags_size_pos,
    const unordered_map<uint32_t, uint32_t>& diags_size_neg,
    vector<pair<uint32_t, char> >& top_diags) {
  /* select the top diags */
  priority_queue<DiagSize> top_diags_queue;
  uint32_t num_top_diags_threshold = num_top_diags + 1;

  /* positive strand */
  for (unordered_map<uint32_t, uint32_t>::const_iterator it = diags_size_pos
      .begin(); it != diags_size_pos.end(); it++) {
    top_diags_queue.push(DiagSize(it->first, it->second, '+'));
    if (top_diags_queue.size() == num_top_diags_threshold) {
      top_diags_queue.pop();
    }
  }

  /* negative strand */
  for (unordered_map<uint32_t, uint32_t>::const_iterator it = diags_size_neg
      .begin(); it != diags_size_neg.end(); it++) {
    top_diags_queue.push(DiagSize(it->first, it->second, '-'));
    if (top_diags_queue.size() == num_top_diags_threshold) {
      top_diags_queue.pop();
    }
  }

  while (!top_diags_queue.empty()) {
    top_diags.push_back(
        make_pair(top_diags_queue.top().diag, top_diags_queue.top().strand));
    top_diags_queue.pop();
  }
}

void Mapping::SingleEndMapping(const string& read) {
  uint32_t read_len = read.size();
  uint32_t num_of_seeds = read_len - HASHLEN + 1;
  string reverse_complement_read = ReverseComplimentStrand(read);
  unordered_map <uint32_t, uint32_t> diags_size_pos;
  unordered_map <uint32_t, uint32_t> diags_size_neg;
  GetDiagsSize(read.c_str(), read_len, diags_size_pos);
  GetDiagsSize(reverse_complement_read.c_str(), read_len, diags_size_neg);

  vector<pair<uint32_t, char> > top_diags;
  GetTopDiags(diags_size_pos, diags_size_neg, top_diags);

  /* select the best match */
  BestMatch best_match;

  int top_diags_size = top_diags.size();
  for (int i = top_diags_size - 1; i >= 0; --i) {
    uint32_t diag_size = 0;
    if (top_diags[i].second == '+') {
      diag_size = diags_size_pos[top_diags[i].first];
    } else {
      diag_size = diags_size_neg[top_diags[i].first];
    }

    uint32_t start_pos = GetGenomeStartPos(top_diags[i].first, read_len);
    if (start_pos + read_len >= genome->all_chroms_len)
      continue;

    if (diag_size == num_of_seeds) {
      /* exact match */
      if (best_match.mismatch > 0) {
        best_match = BestMatch(top_diags[i].first, 1, 0, start_pos,
                               top_diags[i].second);
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
      if (genome->chrom_seqs[k] != read[j]) {
        num_of_mismatch++;
      }
      if (num_of_mismatch > best_match.mismatch)
        break;
      k++;
    }

    if (num_of_mismatch < best_match.mismatch) {
      best_match = BestMatch(top_diags[i].first, 1, num_of_mismatch, start_pos,
                             top_diags[i].second);
    } else if (best_match.mismatch == num_of_mismatch) {
      best_match.times++;
    }
  }
  if (best_match.times > 1) {
    cerr << "[best match to multiple positions, ignore this read]" << endl;
    return;
  }
}
