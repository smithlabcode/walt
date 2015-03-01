#include <utility>
#include <string>
#include <vector>

#include "mapping.hpp"

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

uint32_t LowerBound(uint32_t low, uint32_t high, const char& chr,
                    const uint32_t& cmp_pos,
                    const vector<GenomePosition>& positions,
                    const Genome& genome) {
  uint32_t mid = 0;
  while (low < high) {
    mid = (low + high) / 2;
    char c = genome[positions[mid].chrom_id].sequence[positions[mid].chrom_pos
        + F2SEEDPOSITION[cmp_pos]];
    if (c >= chr) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }
  return low;
}

uint32_t UpperBound(uint32_t low, uint32_t high, const char& chr,
                    const uint32_t& cmp_pos,
                    const vector<GenomePosition>& positions,
                    const Genome& genome) {
  uint32_t mid = 0;
  while (low < high) {
    mid = (low + high + 1) / 2;
    char c = genome[positions[mid].chrom_id].sequence[positions[mid].chrom_pos
        + F2SEEDPOSITION[cmp_pos]];
    if (c <= chr) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }
  return low;
}

void GetRegion(const string& read, const vector<GenomePosition>& positions,
               const Genome& genome, const uint32_t& seed_length,
               pair<uint32_t, uint32_t>& region, TEST_TIME& test_time) {
  test_time.get_region_start_t = clock();
  region.first = 1;
  region.second = 0;

  if (positions.size() == 0)
    return;

  uint32_t l = 0, u = positions.size() - 1;
  for (uint32_t p = F2SEEDWIGTH; p < seed_length; ++p) {
    l = LowerBound(l, u, read[F2SEEDPOSITION[p]], p, positions, genome);
    u = UpperBound(l, u, read[F2SEEDPOSITION[p]], p, positions, genome);
  }
  test_time.get_region_start_sum_time +=
      (clock() - test_time.get_region_start_t);
  if (l > u)
    return;

  region.first = l;
  region.second = u;
}

void SingleEndMapping(const string& orginal_read, const Genome& genome,
                      const HashTable& hash_table, BestMatch& best_match,
                      const uint32_t& seed_length, TEST_TIME& test_time) {
  uint32_t read_len = orginal_read.size();
  if (read_len < HASHLEN)
    return;

  unordered_set<uint64_t> checked_pos;
  string read;
  C2T(orginal_read, read_len, read);
  for (uint32_t seed_i = 0; seed_i < 7; ++seed_i) {
    if (best_match.mismatch == 0)
      break;
    string read_seed = read.substr(seed_i);
    uint32_t hash_value = getHashValue(read_seed.c_str());
    HashTable::const_iterator it = hash_table.find(hash_value);
    if (it == hash_table.end())
      continue;

    pair<uint32_t, uint32_t> region;
    GetRegion(read_seed, it->second, genome, seed_length, region, test_time);
    for (uint32_t j = region.first; j <= region.second; ++j) {
      if (it->second[j].chrom_pos < seed_i)
        continue;
      uint32_t chrom_pos = it->second[j].chrom_pos - seed_i;
      const Chromosome& chrom = genome[it->second[j].chrom_id];
      if (chrom_pos + read_len >= chrom.length)
        continue;

      /* check whether this position already did full check */
      uint64_t check_pos = it->second[j].chrom_id;
      check_pos <<= 32;
      check_pos += chrom_pos;
      if (checked_pos.find(check_pos) != checked_pos.end())
        continue;
      checked_pos.insert(check_pos);

      /* check the position */
      test_time.full_check_time_t = clock();
      uint32_t num_of_mismatch = 0;
      for (uint32_t q = chrom_pos, p = 0;
          p < read_len && num_of_mismatch <= best_match.mismatch; ++q, ++p) {
        if (chrom.sequence[q] != read[p]) {
          num_of_mismatch++;
        }
      }
      test_time.num_of_full_check++;
      if (num_of_mismatch < best_match.mismatch) {
        best_match = BestMatch(it->second[j].chrom_id, chrom_pos, 1,
                               num_of_mismatch);
      } else if (best_match.mismatch == num_of_mismatch) {
        if (best_match.chrom_id != it->second[j].chrom_id
            || best_match.chrom_pos != chrom_pos) {
          best_match.chrom_id = it->second[j].chrom_id;
          best_match.chrom_pos = chrom_pos;
          best_match.times++;
        }
      }
      test_time.full_check_sum_time += (clock() - test_time.full_check_time_t);
    }
  }
}
