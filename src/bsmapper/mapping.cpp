#include <utility>
#include <string>
#include <vector>
#include <algorithm>

#include "mapping.hpp"

using std::count;
using std::copy;

uint32_t MAX(const uint32_t& a, const uint32_t& b) {
  return a > b ? a : b;
}

uint32_t MIN(const uint32_t& a, const uint32_t& b) {
  return a < b ? a : b;
}

string ReverseString(const string& str) {
  uint32_t size = str.size();
  string ret(size, 'N');
  for (uint32_t i = 0; i < size; ++i) {
    ret[i] = str[size - i - 1];
  }

  return ret;
}

string ReverseComplimentString(const string& str) {
  string ret = ReverseString(str);
  for (uint32_t i = 0; i < str.size(); ++i) {
    ret[i] = complimentBase(ret[i]);
  }

  return ret;
}

void ForwardGenomePosition(const uint32_t& genome_pos, const char& strand,
                           const uint32_t& chr_id, const uint32_t& read_length,
                           const Genome& genome, uint32_t& s, uint32_t& e) {
  s = genome_pos - genome.start_index[chr_id];
  s = strand == '+' ? s : genome.length[chr_id] - s - read_length;
  e = s + read_length;
}

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
          p < read_len && num_of_mismatch <= max_mismatches; ++q, ++p) {
        if (genome.sequence[q] != read[p]) {
          num_of_mismatch++;
        }
      }

      if (num_of_mismatch > max_mismatches) {
        continue;
      }
      test_time.num_of_full_check++;
      top_match.Push(CandidatePosition(genome_pos, strand, num_of_mismatch));
      test_time.full_check_sum_time += (clock() - test_time.full_check_time_t);
    }
  }
}

void OutputBestPairedResults(const CandidatePosition& r1,
                             const CandidatePosition& r2, const int& frag_range,
                             const uint32_t& read_length, const Genome& genome,
                             const string& read_name, const string& read_seq1,
                             const string& read_score1, const string& read_seq2,
                             const string& read_score2, ofstream& fout) {

  string read_seq2_rev = ReverseComplimentString(read_seq2);
  string read_scr2_rev = ReverseString(read_score2);

  uint32_t chr_id1 = getChromID(genome.start_index, r1.genome_pos);
  uint32_t chr_id2 = getChromID(genome.start_index, r2.genome_pos);

  uint32_t s1 = 0, s2 = 0, e1 = 0, e2 = 0;
  ForwardGenomePosition(r1.genome_pos, r1.strand, chr_id1, read_length, genome,
                        s1, e1);
  ForwardGenomePosition(r2.genome_pos, r2.strand, chr_id2, read_length, genome,
                        s2, e2);

  uint32_t overlap_s = MAX(s1, s2);
  uint32_t overlap_e = MIN(e1, e2);

  uint32_t one_l = r1.strand == '+' ? s1 : MAX(overlap_e, s1);
  uint32_t one_r = r1.strand == '+' ? MIN(overlap_s, e1) : e1;

  uint32_t two_l = r1.strand == '+' ? MAX(overlap_e, s2) : s2;
  uint32_t two_r = r1.strand == '+' ? e2 : MIN(overlap_s, e2);

  int len = r1.strand == '+' ? (two_r - one_l) : (one_r - two_l);

  string seq(len, 'N');
  string scr(len, 'B');
  if (len > 0 && len <= frag_range) {
    // lim_one: offset in merged sequence where overlap starts
    uint32_t lim_one = one_r - one_l;
    copy(read_seq1.begin(), read_seq1.begin() + lim_one, seq.begin());
    copy(read_score1.begin(), read_score1.begin() + lim_one, scr.begin());

    uint32_t lim_two = two_r - two_l;
    copy(read_seq2_rev.end() - lim_two, read_seq2_rev.end(),
         seq.end() - lim_two);
    copy(read_scr2_rev.end() - lim_two, read_scr2_rev.end(),
         scr.end() - lim_two);

    // deal with overlapping part
    if (overlap_s < overlap_e) {
      uint32_t one_bads = count(read_seq1.begin(), read_seq1.end(), 'N');
      int info_one = read_length - (one_bads + r1.mismatch);

      uint32_t two_bads = count(read_seq2_rev.begin(), read_seq2_rev.end(),
                                'N');
      int info_two = read_length - (two_bads + r2.mismatch);

      // use the mate with the most info to fill in the overlap
      if (info_one >= info_two) {
        uint32_t a = r1.strand == '+' ? (overlap_s - s1) : (e1 - overlap_e);
        uint32_t b = r1.strand == '+' ? (overlap_e - s1) : (e1 - overlap_s);
        copy(read_seq1.begin() + a, read_seq1.begin() + b,
             seq.begin() + lim_one);
        copy(read_score1.begin() + a, read_score1.begin() + b,
             scr.begin() + lim_one);
      } else {
        uint32_t a = r1.strand == '+' ? (overlap_s - s2) : (e2 - overlap_e);
        uint32_t b = r1.strand == '+' ? (overlap_e - s2) : (e2 - overlap_s);
        copy(read_seq2_rev.begin() + a, read_seq2_rev.begin() + b,
             seq.begin() + lim_one);
        copy(read_scr2_rev.begin() + a, read_scr2_rev.begin() + b,
             scr.begin() + lim_one);
      }
    }
  }

  uint32_t start_pos = r1.strand == '+' ? s1 : s2;
  fout << genome.name[chr_id1] << "\t" << start_pos << "\t" << start_pos + len
      << "\t" << "FRAG:" << read_name << "\t" << r1.mismatch + r2.mismatch
      << "\t" << r1.strand << "\t" << seq << "\t" << scr << endl;
}

void OutputBestSingleResults(const vector<CandidatePosition>& ranked_results,
                             const int ranked_results_size,
                             const Genome& genome, const uint32_t& read_length,
                             const string& read_name, const string& read_seq,
                             const string& read_score,
                             const uint32_t& max_mismatches, ofstream& fout) {
  BestMatch best_match(0, 0, '+', max_mismatches);
  for (int i = ranked_results_size - 1; i >= 0; --i) {
    const CandidatePosition& r = ranked_results[i];
    if (r.mismatch < best_match.mismatch) {
      best_match = BestMatch(r.genome_pos, 1, r.strand, r.mismatch);
    } else if (r.mismatch == best_match.mismatch
        && best_match.genome_pos != r.genome_pos) {
      best_match.genome_pos = r.genome_pos;
      best_match.strand = r.strand;
      best_match.times++;
    } else {
      break;
    }
  }

  if (best_match.times == 1) {
    uint32_t chr_id = getChromID(genome.start_index, best_match.genome_pos);
    uint32_t start_pos = 0, end_pos = 0;
    ForwardGenomePosition(best_match.genome_pos, best_match.strand, chr_id,
                          read_length, genome, start_pos, end_pos);

    fout << genome.name[chr_id] << "\t" << start_pos << "\t" << end_pos << "\t"
        << read_name << "\t" << best_match.mismatch << "\t" << best_match.strand
        << "\t" << read_seq << "\t" << read_score << endl;
  }
}

int GetFragmentLength(const CandidatePosition& r1, const CandidatePosition& r2,
                      const uint32_t& frag_range, const uint32_t& read_length,
                      const Genome& genome, const uint32_t& chr_id1,
                      const uint32_t& chr_id2) {
  uint32_t s1 = 0, s2 = 0, e1 = 0, e2 = 0;
  ForwardGenomePosition(r1.genome_pos, r1.strand, chr_id1, read_length, genome,
                        s1, e1);
  ForwardGenomePosition(r2.genome_pos, r2.strand, chr_id2, read_length, genome,
                        s2, e2);

  return r1.strand == '+' ? (e2 - s1) : (e1 - s2);
}

void MergePairedEndResults(
    const vector<vector<CandidatePosition> >& ranked_results,
    const vector<int>& ranked_results_size, const uint32_t& max_mismatches,
    const uint32_t& read_length, const int& frag_range, const Genome& genome,
    const string& read_name, const string& read_seq1, const string& read_score1,
    const string& read_seq2, const string& read_score2, ofstream& fout) {
#ifdef DEBUG
  for (int i = ranked_results_size[0] - 1; i >= 0; --i) {
    const CandidatePosition& r = ranked_results[0][i];
    cerr << "LL " << i << ": " << r.genome_pos << " " << r.strand << " "
    << r.mismatch << endl;
  }
  for (int i = ranked_results_size[1] - 1; i >= 0; --i) {
    const CandidatePosition& r = ranked_results[1][i];
    cerr << "UU " << i << ": " << r.genome_pos << " " << r.strand << " "
    << r.mismatch << endl;
  }
#endif

  pair<int, int> best_pair(-1, -1);
  uint32_t min_num_of_mismatch = max_mismatches;
  uint64_t best_pos = 0;
  uint32_t best_times = 0;
  for (int i = ranked_results_size[0] - 1; i >= 0; --i) {
    for (int j = ranked_results_size[1] - 1; j >= 0; --j) {
      const CandidatePosition& r1 = ranked_results[0][i];
      const CandidatePosition& r2 = ranked_results[1][j];
      if (r1.strand == r2.strand)
        continue;

      uint32_t num_of_mismatch = r1.mismatch + r2.mismatch;
      if (num_of_mismatch > min_num_of_mismatch)
        break;

      uint32_t chr_id1 = getChromID(genome.start_index, r1.genome_pos);
      uint32_t chr_id2 = getChromID(genome.start_index, r2.genome_pos);
      if (chr_id1 != chr_id2)
        continue;

      int frag_size = GetFragmentLength(r1, r2, frag_range, read_length, genome,
                                        chr_id1, chr_id2);
      if (frag_size <= 0 || frag_size > frag_range)
        continue;

      uint64_t cur_pos = r1.genome_pos;
      cur_pos <<= 32;
      cur_pos += r2.genome_pos;
      if (num_of_mismatch < min_num_of_mismatch) {
        best_pair = make_pair(i, j);
        best_times = 1;
        min_num_of_mismatch = num_of_mismatch;
        best_pos = cur_pos;
      } else if (num_of_mismatch == min_num_of_mismatch
          && cur_pos != best_pos) {
        best_pair = make_pair(i, j);
        best_times++;
      }
    }
  }

  if (best_times == 1) {
    OutputBestPairedResults(ranked_results[0][best_pair.first],
                            ranked_results[1][best_pair.second], frag_range,
                            read_length, genome, read_name, read_seq1,
                            read_score1, read_seq2, read_score2, fout);
  } else {
    OutputBestSingleResults(ranked_results[0], ranked_results_size[0], genome,
                            read_length, read_name, read_seq1, read_score1,
                            max_mismatches, fout);
    OutputBestSingleResults(ranked_results[1], ranked_results_size[1], genome,
                            read_length, read_name, read_seq2, read_score2,
                            max_mismatches, fout);
  }
}
