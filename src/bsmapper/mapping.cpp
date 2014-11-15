#include "mapping.hpp"

#include <limits>
#include <queue>
#include <fstream>

#include <tr1/unordered_map>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::tr1::unordered_map;
using std::priority_queue;

void LoadReadsFromFastqFile(const string &filename,
                            const uint64_t read_start_idx,
                            const uint64_t n_reads_to_process,
                            vector<string>& read_names,
                            vector<string>& read_seqs) {
  if (n_reads_to_process != std::numeric_limits<uint64_t>::max()) {
    cerr << "[LOADING READS FROM " << read_start_idx << " TO "
         << n_reads_to_process + read_start_idx << "]" << endl;
  } else {
    cerr << "[LOADING READS FROM " << read_start_idx << " TO LAST ONE]" << endl;
  }
  read_names.clear();
  read_seqs.clear();
  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + filename);

  uint64_t line_count = 0;
  const uint64_t lim1 = read_start_idx * 4;
  while (line_count < lim1) {
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    ++line_count;
  }

  const uint64_t lim2 =
      (n_reads_to_process != std::numeric_limits<size_t>::max()) ?
          (read_start_idx + n_reads_to_process) * 4 :
          std::numeric_limits<size_t>::max();

  string line;
  while (line_count < lim2 && getline(in, line)) {
    if (isFastqSequenceLine(line_count)) {
      read_seqs.push_back(line);
    } else if (isFastqNameLine(line_count)) {
      read_names.push_back(line);
    }
    ++line_count;
  }
}

inline uint32_t getDiag(const uint32_t& genome_pos, const uint32_t& seed_pos,
                        const uint32_t& read_len) {
  return genome_pos + read_len - seed_pos;
}

inline uint32_t getGenomeStartPos(const uint32_t& diag,
                                  const uint32_t& read_len) {
  if (read_len > diag)
    return 0;
  return diag - read_len;
}

void mapping(const Genome* genome, const HashTable* hash_table,
             const char* read, const int& num_top_diags) {
  /* count how many matched seeds in each diag */
  unordered_map<uint32_t, uint32_t> diags_size;
  uint32_t read_len = strlen(read), hash_value = 0;
  uint32_t num_of_seeds = read_len - HASHLEN + 1;
  for (uint32_t i = 0; i < read_len - HASHLEN + 1; ++i) {
    hash_value = getHashValue(&(read[i]));
    pair<HashTable::const_iterator, HashTable::const_iterator> range = hash_table->equal_range(hash_value);
    for (HashTable::const_iterator val = range.first; val != range.second; ++val) {
      diags_size[getDiag(val->second, i, read_len)]++;
    }
  }

  /* select the top diags */
  priority_queue<DiagSize> top_diags_queue;
  uint32_t num_top_diags_threshold = num_top_diags + 1;
  for (unordered_map<uint32_t, uint32_t>::const_iterator it =
      diags_size.begin(); it != diags_size.end(); it++) {
    top_diags_queue.push(DiagSize(it->first, it->second));
    if (top_diags_queue.size() == num_top_diags_threshold) {
      top_diags_queue.pop();
    }
  }
  vector<uint32_t> top_diags;
  while (!top_diags_queue.empty()) {
    top_diags.push_back(top_diags_queue.top().diag);
    top_diags_queue.pop();
  }

  /* select the best match */
  uint32_t best_match_diag;
  uint32_t best_match_diag_times = 0;
  uint32_t best_match_mismatch = std::numeric_limits<uint32_t>::max();

  int top_diags_size = top_diags.size();
  for (int i = top_diags_size - 1; i >= 0; --i) {
    if (diags_size[top_diags[i]] == num_of_seeds) {
      /* exact match */
      if (best_match_mismatch > 0) {
        best_match_diag = top_diags[i];
        best_match_mismatch = 0;
        best_match_diag_times = 1;
        continue;
      }

      if (best_match_mismatch == 0) {
        cerr << "[map to multiple positions, ignore this read]" << endl;
        return;
      }
    }

    if (best_match_mismatch == 0)
      continue;

    /* mismatch */
    uint32_t start_pos = getGenomeStartPos(top_diags[i], read_len);
    if (start_pos + read_len >= genome->all_chroms_len)
      continue;

    uint32_t num_of_mismatch = 0;
    for (uint32_t j = 0; j < read_len; ++j) {
      if (genome->chrom_seqs[start_pos] != read[j]) {
        num_of_mismatch++;
      }
      if (num_of_mismatch > best_match_mismatch)
        break;
      start_pos++;
    }

    if (num_of_mismatch < best_match_mismatch) {
      best_match_diag = top_diags[i];
      best_match_mismatch = num_of_mismatch;
      best_match_diag_times = 1;
    } else if (best_match_mismatch == num_of_mismatch) {
      best_match_diag_times++;
    }
  }
  if (best_match_diag_times > 1) {
    cerr << "[map to multiple positions, ignore this read]" << endl;
    return;
  }
}
