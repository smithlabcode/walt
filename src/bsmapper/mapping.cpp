#include <utility>
#include <string>
#include <vector>

#include "mapping.hpp"

/* for bisulfite sequence mapping, Cs are transfered to Ts*/
void C2T(const string& orginal_read, const uint32_t& read_len, string& read) {
  for (uint32_t i = 0; i < read_len; ++i) {
    if ('N' == orginal_read[i]) {
      read += 'T';  // in rmapbs N set to 3.
    } else if ('C' == orginal_read[i]) {
      read += 'T';
    } else {
      read += orginal_read[i];
    }
  }
}

void A2G(const string& orginal_read, const uint32_t& read_len, string& read) {
  for (uint32_t i = 0; i < read_len; ++i) {
    if ('N' == orginal_read[i]) {
      read += 'G';  // in rmapbs N set to 3.
    } else if ('A' == orginal_read[i]) {
      read += 'G';
    } else {
      read += orginal_read[i];
    }
  }
}

uint32_t LowerBound(uint32_t low, uint32_t high, const char& chr,
                    const uint32_t& cmp_pos, const Genome& genome,
                    const HashTable& hash_table) {
  uint32_t mid = 0;
  while (low < high) {
    mid = low + (high - low) / 2;
    char c = genome.sequence[hash_table.index[mid] + cmp_pos];
    if (c >= chr) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }
  return low;
}

uint32_t UpperBound(uint32_t low, uint32_t high, const char& chr,
                    const uint32_t& cmp_pos, const Genome& genome,
                    const HashTable& hash_table) {
  uint32_t mid = 0;
  while (low < high) {
    mid = low + (high - low + 1) / 2;
    char c = genome.sequence[hash_table.index[mid] + cmp_pos];
    if (c <= chr) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }
  return low;
}

void GetRegion(const string& read, const Genome& genome,
               const HashTable& hash_table, const uint32_t& seed_length,
               pair<uint32_t, uint32_t>& region, TEST_TIME& test_time) {
  test_time.get_region_start_t = clock();

  uint32_t l = region.first, u = region.second - 1;
  for (uint32_t p = F2SEEDWIGTH; p < seed_length; ++p) {
    uint32_t care_pos = F2SEEDPOSITION[p];
    l = LowerBound(l, u, read[care_pos], care_pos, genome, hash_table);
    u = UpperBound(l, u, read[care_pos], care_pos, genome, hash_table);
  }
  test_time.get_region_start_sum_time +=
      (clock() - test_time.get_region_start_t);
  if (l > u) {
    region.first = 1;
    region.second = 0;
    return;
  }

  region.first = l;
  region.second = u;
}

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      const HashTable& hash_table, BestMatch& best_match,
                      const uint32_t& seed_length, const char& strand,
                      TEST_TIME& test_time, const bool& AG_WILDCARD) {
  uint32_t read_len = orginal_read.size();

  string read;
  if (AG_WILDCARD) {
    A2G(orginal_read, read_len, read);
  } else {
    C2T(orginal_read, read_len, read);
  }

  for (uint32_t seed_i = 0; seed_i < SEEPATTERNLEN; ++seed_i) {
    if (best_match.mismatch == 0)
      break;
    string read_seed = read.substr(seed_i);
    uint32_t hash_value = getHashValue(read_seed.c_str());
    pair<uint32_t, uint32_t> region;
    region.first = hash_table.counter[hash_value];
    region.second = hash_table.counter[hash_value + 1];

    if (region.first == region.second)
      continue;

    GetRegion(read_seed, genome, hash_table, seed_length, region, test_time);
    for (uint32_t j = region.first; j <= region.second; ++j) {
      uint32_t genome_pos = hash_table.index[j];
      uint32_t chr_id = getChromID(genome.start_index, genome_pos);
      if (genome_pos - genome.start_index[chr_id] < seed_i)
        continue;
      genome_pos = genome_pos - seed_i;
      if (genome_pos + read_len >= genome.start_index[chr_id + 1])
        continue;

      /* check the position */
      test_time.full_check_time_t = clock();
      uint32_t num_of_mismatch = 0;
      for (uint32_t q = genome_pos, p = 0;
          p < read_len && num_of_mismatch <= best_match.mismatch; ++q, ++p) {
        if (genome.sequence[q] != read[p]) {
          num_of_mismatch++;
        }
      }
      test_time.num_of_full_check++;
      if (num_of_mismatch < best_match.mismatch) {
        best_match = BestMatch(genome_pos, 1, strand, num_of_mismatch);
      } else if (best_match.mismatch == num_of_mismatch
          && best_match.genome_pos != genome_pos) {
        best_match.genome_pos = genome_pos;
        best_match.strand = strand;
        best_match.times++;
      }
      test_time.full_check_sum_time += (clock() - test_time.full_check_time_t);
    }
  }
}

void PairEndMapping(const string& orginal_read, const Genome& genome,
                    const HashTable& hash_table, TopCandidates& top_match,
                    const uint32_t& max_mismatches, const uint32_t& seed_length,
                    const char& strand, const bool& AG_WILDCARD) {
  uint32_t read_len = orginal_read.size();

  string read;
  if (AG_WILDCARD) {
    A2G(orginal_read, read_len, read);
  } else {
    C2T(orginal_read, read_len, read);
  }

  for (uint32_t seed_i = 0; seed_i < SEEPATTERNLEN; ++seed_i) {
    if (!top_match.Empty() && top_match.Top().mismatch == 0)
      break;

    string read_seed = read.substr(seed_i);
    uint32_t hash_value = getHashValue(read_seed.c_str());
    pair<uint32_t, uint32_t> region;
    region.first = hash_table.counter[hash_value];
    region.second = hash_table.counter[hash_value + 1];

    if (region.first == region.second)
      continue;
    TEST_TIME test_time;
    GetRegion(read_seed, genome, hash_table, seed_length, region, test_time);
    uint32_t top_mismatch;
    if (!top_match.Empty()) {
      top_mismatch = top_match.Top().mismatch;
    } else {
      top_mismatch = max_mismatches;
    }
    for (uint32_t j = region.first; j <= region.second; ++j) {
      uint32_t genome_pos = hash_table.index[j];
      uint32_t chr_id = getChromID(genome.start_index, genome_pos);
      if (genome_pos - genome.start_index[chr_id] < seed_i)
        continue;
      genome_pos = genome_pos - seed_i;
      if (genome_pos + read_len >= genome.start_index[chr_id + 1])
        continue;

      /* check the position */
      test_time.full_check_time_t = clock();
      uint32_t num_of_mismatch = 0;
      for (uint32_t q = genome_pos, p = 0;
          p < read_len && num_of_mismatch <= top_mismatch; ++q, ++p) {
        if (genome.sequence[q] != read[p]) {
          num_of_mismatch++;
        }
      }

      test_time.num_of_full_check++;
      top_match.Push(CandidatePosition(genome_pos, strand, num_of_mismatch));
      test_time.full_check_sum_time += (clock() - test_time.full_check_time_t);
    }
  }
}

void MergePairedEndResults(
    const vector<vector<CandidatePosition> >& ranked_results,
    const vector<int>& ranked_results_size, BestMatch& best_match,
    const uint32_t& L, const uint32_t& U, const Genome& genome) {
  for (int i = ranked_results_size[0] - 1; i >= 0; --i) {
    for (int j = ranked_results_size[1]; j >= 0; --j) {
      const CandidatePosition& r1 = ranked_results[0][i];
      const CandidatePosition& r2 = ranked_results[1][j];
      uint32_t num_of_mismatch = r1.mismatch + r2.mismatch;
      if (num_of_mismatch > best_match.mismatch)
        break;

      uint32_t dis =
          r1.genome_pos >= r2.genome_pos ?
              r1.genome_pos - r2.genome_pos : r2.genome_pos - r1.genome_pos;
      if (r1.strand != r2.strand && dis >= L && dis <= U) {
        uint32_t chr_id1 = getChromID(genome.start_index, r1.genome_pos);
        uint32_t chr_id2 = getChromID(genome.start_index, r2.genome_pos);
        if (chr_id1 != chr_id2)
          continue;

        uint32_t genome_pos = r1.strand == '+' ? r1.genome_pos : r2.genome_pos;
        if (num_of_mismatch < best_match.mismatch) {
          best_match = BestMatch(genome_pos, 1, '+', num_of_mismatch);
        } else if (best_match.mismatch == num_of_mismatch
            && best_match.genome_pos != genome_pos) {
          best_match.genome_pos = genome_pos;
          best_match.strand = '+';
          best_match.times++;
        }
      }
    }
  }
}
