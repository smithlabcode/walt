#include "paired.hpp"

#include <algorithm>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

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

/* find the genome position for the in the forward strand */
void ForwardChromPosition(const uint32_t& genome_pos, const char& strand,
                          const uint32_t& chr_id, const uint32_t& read_len,
                          const Genome& genome, uint32_t& s, uint32_t& e) {
  s = genome_pos - genome.start_index[chr_id];
  s = strand == '+' ? s : genome.length[chr_id] - s - read_len;
  e = s + read_len;
}

void PairEndMapping(const string& org_read, const Genome& genome,
                    const HashTable& hash_table, const char& strand,
                    const bool& AG_WILDCARD, const uint32_t& max_mismatches,
                    const uint32_t& seed_len, TopCandidates& top_match) {
  uint32_t read_len = org_read.size();

  string read;
  if (AG_WILDCARD) {
    A2G(org_read, read_len, read);
  } else {
    C2T(org_read, read_len, read);
  }

  for (uint32_t seed_i = 0; seed_i < SEEPATTERNLEN; ++seed_i) {
    string read_seed = read.substr(seed_i);
    uint32_t hash_value = getHashValue(read_seed.c_str());
    pair < uint32_t, uint32_t > region;
    region.first = hash_table.counter[hash_value];
    region.second = hash_table.counter[hash_value + 1];

    if (region.first == region.second)
      continue;

    IndexRegion(read_seed, genome, hash_table, seed_len, region);
    for (uint32_t j = region.first; j <= region.second; ++j) {
      uint32_t genome_pos = hash_table.index[j];
      uint32_t chr_id = getChromID(genome.start_index, genome_pos);
      if (genome_pos - genome.start_index[chr_id] < seed_i)
        continue;
      genome_pos = genome_pos - seed_i;
      if (genome_pos + read_len >= genome.start_index[chr_id + 1])
        continue;

      /* check the position */
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
      top_match.Push(CandidatePosition(genome_pos, strand, num_of_mismatch));
    }
  }
}

void OutputBestPairedResults(const CandidatePosition& r1,
                             const CandidatePosition& r2, const int& frag_range,
                             const uint32_t& read_len, const Genome& genome,
                             const string& read_name, const string& read_seq1,
                             const string& read_score1, const string& read_seq2,
                             const string& read_score2, ofstream& fout) {

  string read_seq2_rev = ReverseComplimentString(read_seq2);
  string read_scr2_rev = ReverseString(read_score2);

  uint32_t chr_id1 = getChromID(genome.start_index, r1.genome_pos);
  uint32_t chr_id2 = getChromID(genome.start_index, r2.genome_pos);

  uint32_t s1 = 0, s2 = 0, e1 = 0, e2 = 0;
  ForwardChromPosition(r1.genome_pos, r1.strand, chr_id1, read_len, genome, s1,
                       e1);
  ForwardChromPosition(r2.genome_pos, r2.strand, chr_id2, read_len, genome, s2,
                       e2);

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
      int info_one = read_len - (one_bads + r1.mismatch);

      uint32_t two_bads = count(read_seq2_rev.begin(), read_seq2_rev.end(),
                                'N');
      int info_two = read_len - (two_bads + r2.mismatch);

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
                             const Genome& genome, const uint32_t& read_len,
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
    ForwardChromPosition(best_match.genome_pos, best_match.strand, chr_id,
                         read_len, genome, start_pos, end_pos);

    fout << genome.name[chr_id] << "\t" << start_pos << "\t" << end_pos << "\t"
        << read_name << "\t" << best_match.mismatch << "\t" << best_match.strand
        << "\t" << read_seq << "\t" << read_score << endl;
  }
}

int GetFragmentLength(const CandidatePosition& r1, const CandidatePosition& r2,
                      const uint32_t& frag_range, const uint32_t& read_len,
                      const Genome& genome, const uint32_t& chr_id1,
                      const uint32_t& chr_id2) {
  uint32_t s1 = 0, s2 = 0, e1 = 0, e2 = 0;
  ForwardChromPosition(r1.genome_pos, r1.strand, chr_id1, read_len, genome, s1,
                       e1);
  ForwardChromPosition(r2.genome_pos, r2.strand, chr_id2, read_len, genome, s2,
                       e2);

  return r1.strand == '+' ? (e2 - s1) : (e1 - s2);
}

/* merge the mapping results from paired reads */
void MergePairedEndResults(
    const Genome& genome, const string& read_name, const string& read_seq1,
    const string& read_score1, const string& read_seq2,
    const string& read_score2,
    const vector<vector<CandidatePosition> >& ranked_results,
    const vector<int>& ranked_results_size, const uint32_t& read_len,
    const int& frag_range, const uint32_t& max_mismatches, ofstream& fout) {
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

      int frag_size = GetFragmentLength(r1, r2, frag_range, read_len, genome,
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
                            read_len, genome, read_name, read_seq1, read_score1,
                            read_seq2, read_score2, fout);
  } else {
    OutputBestSingleResults(ranked_results[0], ranked_results_size[0], genome,
                            read_len, read_name, read_seq1, read_score1,
                            max_mismatches, fout);
    OutputBestSingleResults(ranked_results[1], ranked_results_size[1], genome,
                            read_len, read_name, read_seq2, read_score2,
                            max_mismatches, fout);
  }
}

void ProcessPairedEndReads(const string& index_file,
                           const uint32_t& n_reads_to_process,
                           const string& reads_file_p1,
                           const string& reads_file_p2,
                           const string& output_file,
                           const uint32_t& max_mismatches,
                           const uint32_t& read_len, const uint32_t& seed_len,
                           const uint32_t& top_k, const int& frag_range) {
  // LOAD THE INDEX HEAD INFO
  Genome genome;
  HashTable hash_table;

  uint32_t size_of_index;
  ReadIndexHeadInfo(index_file, genome, size_of_index);
  genome.sequence.resize(genome.length_of_genome);
  hash_table.counter.resize(power(4, F2SEEDWIGTH) + 1);
  hash_table.index.resize(size_of_index);

  vector<vector<string> > index_names(2, vector < string > (2));
  index_names[0][0] = index_file + "_CT00";
  index_names[0][1] = index_file + "_CT01";
  index_names[1][0] = index_file + "_AG10";
  index_names[1][1] = index_file + "_AG11";

  vector<vector<string> > read_names(2, vector<string>(n_reads_to_process));
  vector<vector<string> > read_seqs(2, vector<string>(n_reads_to_process));
  vector<vector<string> > read_scores(2, vector<string>(n_reads_to_process));

  vector<vector<TopCandidates> > top_results(2,
                vector<TopCandidates>(n_reads_to_process));

  vector<BestMatch> map_results(n_reads_to_process);
  vector<int> ranked_results_size(2);
  vector<vector<CandidatePosition> > ranked_results(2,
                vector<CandidatePosition>(top_k));

  ifstream fin[2];
  fin[0].open(reads_file_p1.c_str());
  fin[1].open(reads_file_p2.c_str());
  if (!fin[0]) {
    throw SMITHLABException("cannot open input file " + reads_file_p1);
  }
  if (!fin[1]) {
    throw SMITHLABException("cannot open input file " + reads_file_p1);
  }

  clock_t start_t;
  uint64_t sum_t = 0;

  ofstream fout(output_file.c_str());
  uint32_t num_of_reads[2];
  bool AG_WILDCARD = true;
  for (uint64_t i = 0;; i += n_reads_to_process) {
    for (uint32_t pi = 0; pi < 2; ++pi) {  // paired end reads _1 and _2
      AG_WILDCARD = pi == 1 ? true : false;

      LoadReadsFromFastqFile(fin[pi], i, n_reads_to_process, num_of_reads[pi],
                             read_names[pi], read_seqs[pi], read_scores[pi]);
      if (num_of_reads[pi] == 0)
        break;
      if (pi == 0) {
        cerr << "[START MAPPING READS FROM " << i << " TO "
            << num_of_reads[pi] + i << "]" << endl;
      }

      //Initialize the paired results
      for (uint32_t j = 0; j < num_of_reads[pi]; ++j) {
        top_results[pi][j].Clear();
        top_results[pi][j].SetSize(top_k);
      }

      for (uint32_t fi = 0; fi < 2; ++fi) {
        TIME_INFO(ReadIndex(index_names[pi][fi], genome, hash_table),
                  "LOAD INDEX");
        for (uint32_t j = 0; j < num_of_reads[pi]; ++j) {
          char strand = fi == 0 ? '+' : '-';
          start_t = clock();
          PairEndMapping(read_seqs[pi][j], genome, hash_table, strand,
                         AG_WILDCARD, max_mismatches, seed_len,
                         top_results[pi][j]);
          sum_t += clock() - start_t;
        }
        fprintf(stderr, "[%.3lf SECONDS MAPPING TIME PASSED]\n",
                static_cast<double>(sum_t / CLOCKS_PER_SEC));
      }
    }

    ///////////////////////////////////////////////////////////
    //Merge Paired-end results
    for (uint32_t j = 0; j < num_of_reads[0]; ++j) {
      DEBUG_INFO(read_names[0][j], "\n");
      for (uint32_t pi = 0; pi < 2; ++pi) {
        ranked_results_size[pi] = 0;
        while (!top_results[pi][j].candidates.empty()) {
          ranked_results[pi][ranked_results_size[pi]++] =
              top_results[pi][j].Top();
          top_results[pi][j].Pop();
        }
      }

      MergePairedEndResults(genome, read_names[0][j], read_seqs[0][j],
                            read_scores[0][j], read_seqs[1][j],
                            read_scores[1][j], ranked_results,
                            ranked_results_size, read_len, frag_range,
                            max_mismatches, fout);
    }

    if (num_of_reads[0] < n_reads_to_process)
      break;
  }

  fin[0].close();
  fin[1].close();
  fout.close();

  fprintf(stderr, "[MAPPING TAKES %.3lf SECONDS]\n",
          static_cast<double>(sum_t / CLOCKS_PER_SEC));
}
