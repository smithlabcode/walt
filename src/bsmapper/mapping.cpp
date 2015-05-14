#include "mapping.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

void LoadReadsFromFastqFile(FILE * fin, const uint32_t read_start_idx,
                            const uint32_t n_reads_to_process,
                            uint32_t& num_of_reads, vector<string>& read_names,
                            vector<string>& read_seqs,
                            vector<string>& read_scores) {
  char cline[MAX_LINE_LENGTH];
  string line;
  int line_code = 0;
  uint32_t line_count = 0;
  num_of_reads = 0;
  uint32_t lim = n_reads_to_process * 4;
  while (line_count < lim && fgets(cline, MAX_LINE_LENGTH, fin)) {
    cline[strlen(cline) - 1] = 0;
    line = cline;
    switch (line_code) {
      case 0: {
        size_t space_pos = line.find_first_of(' ');
        if (space_pos == string::npos) {
          read_names[num_of_reads] = line.substr(1);
        } else {
          read_names[num_of_reads] = line.substr(1, space_pos - 1);
        }
        break;
      }
      case 1: {
        read_seqs[num_of_reads] = line;
        break;
      }
      case 2: {
        break;
      }
      case 3: {
        read_scores[num_of_reads] = line;
        num_of_reads++;
        break;
      }
    }
    ++line_count;
    ++line_code;
    if (line_code == 4) {
      line_code = 0;
    }
  }
}

void C2T(const string& org_read, const uint32_t& read_len, string& read) {
  for (uint32_t i = 0; i < read_len; ++i) {
    if ('N' == org_read[i]) {
      read += 'T';
    } else if ('C' == org_read[i]) {
      read += 'T';
    } else {
      read += org_read[i];
    }
  }
}

void A2G(const string& org_read, const uint32_t& read_len, string& read) {
  for (uint32_t i = 0; i < read_len; ++i) {
    if ('N' == org_read[i]) {
      read += 'G';
    } else if ('A' == org_read[i]) {
      read += 'G';
    } else {
      read += org_read[i];
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

void IndexRegion(const string& read, const Genome& genome,
                 const HashTable& hash_table, const uint32_t& seed_len,
                 pair<uint32_t, uint32_t>& region) {
  uint32_t l = region.first, u = region.second - 1;
  for (uint32_t p = F2SEEDWIGTH; p < seed_len; ++p) {
    uint32_t care_pos = F2SEEDPOSITION[p];
    l = LowerBound(l, u, read[care_pos], care_pos, genome, hash_table);
    u = UpperBound(l, u, read[care_pos], care_pos, genome, hash_table);
  }

  if (l > u) {
    region.first = 1;
    region.second = 0;
    return;
  }

  region.first = l;
  region.second = u;
}

void SingleEndMapping(const string& org_read, const Genome& genome,
                      const HashTable& hash_table, const char& strand,
                      const bool& AG_WILDCARD, const uint32_t& seed_len,
                      BestMatch& best_match) {
  uint32_t read_len = org_read.size();

  string read;
  if (AG_WILDCARD) {
    A2G(org_read, read_len, read);
  } else {
    C2T(org_read, read_len, read);
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
          p < read_len && num_of_mismatch <= best_match.mismatch; ++q, ++p) {
        if (genome.sequence[q] != read[p]) {
          num_of_mismatch++;
        }
      }

      if (num_of_mismatch < best_match.mismatch) {
        best_match = BestMatch(genome_pos, 1, strand, num_of_mismatch);
      } else if (best_match.mismatch == num_of_mismatch
          && best_match.genome_pos != genome_pos) {
        best_match.genome_pos = genome_pos;
        best_match.strand = strand;
        best_match.times++;
      }
    }
  }
}

void OutputSingleEndResults(FILE * fout, const vector<BestMatch>& map_results,
                            const uint32_t& num_of_reads,
                            const vector<string>& read_names,
                            const vector<string>& read_seqs,
                            const vector<string>& read_scores,
                            const Genome& genome) {
  for (uint32_t j = 0; j < num_of_reads; ++j) {
    if (map_results[j].times == 0 || map_results[j].times > 1)
      continue;

    uint32_t chr_id = getChromID(genome.start_index, map_results[j].genome_pos);
    uint32_t start_pos = map_results[j].genome_pos - genome.start_index[chr_id];
    if ('-' == map_results[j].strand) {
      start_pos = genome.length[chr_id] - start_pos - read_seqs[j].size();
    }
    uint32_t end_pos = start_pos + read_seqs[j].size();

    fprintf(fout, "%s\t%u\t%u\t%s\t%u\t%c\t%s\t%s\n",
            genome.name[chr_id].c_str(), start_pos, end_pos,
            read_names[j].c_str(), map_results[j].mismatch,
            map_results[j].strand, read_seqs[j].c_str(),
            read_scores[j].c_str());
  }
}

void ProcessSingledEndReads(const string& index_file,
                            const uint32_t& n_reads_to_process,
                            const string& reads_file_s,
                            const string& output_file,
                            const uint32_t& max_mismatches,
                            const uint32_t& read_len, const uint32_t& seed_len,
                            const bool& AG_WILDCARD) {
  // LOAD THE INDEX HEAD INFO
  Genome genome;
  HashTable hash_table;

  uint32_t size_of_index;
  ReadIndexHeadInfo(index_file, genome, size_of_index);
  genome.sequence.resize(genome.length_of_genome);
  hash_table.counter.resize(power(4, F2SEEDWIGTH) + 1);
  hash_table.index.resize(size_of_index);

  vector<string> index_names;
  if (!AG_WILDCARD) {
    index_names.push_back(index_file + "_CT00");
    index_names.push_back(index_file + "_CT01");
  } else {
    index_names.push_back(index_file + "_AG10");
    index_names.push_back(index_file + "_AG11");
  }

  vector<string> read_names(n_reads_to_process);
  vector<string> read_seqs(n_reads_to_process);
  vector<string> read_scores(n_reads_to_process);

  vector<BestMatch> map_results(n_reads_to_process);
  FILE * fin = fopen(reads_file_s.c_str(), "r");
  if (!fin) {
    throw SMITHLABException("cannot open input file " + reads_file_s);
  }
  clock_t start_t;
  uint64_t sum_t = 0;

  FILE * fout = fopen(output_file.c_str(), "w");
  uint32_t num_of_reads;
  for (uint32_t i = 0;; i += n_reads_to_process) {
    LoadReadsFromFastqFile(fin, i, n_reads_to_process, num_of_reads, read_names,
                           read_seqs, read_scores);
    if (num_of_reads == 0)
      break;

    // Initialize the results
    BestMatch best_match(0, 0, '+', max_mismatches);
    for (uint32_t j = 0; j < num_of_reads; ++j) {
      map_results[j] = best_match;
    }

    fprintf(stderr, "[START MAPPING READS FROM %u TO %u]\n", i,
            num_of_reads + i);
    for (uint32_t fi = 0; fi < 2; ++fi) {
      TIME_INFO(ReadIndex(index_names[fi], genome, hash_table), "LOAD INDEX");
      for (uint32_t j = 0; j < num_of_reads; ++j) {
        DEBUG_INFO(read_names[j], "\n");
        char strand = fi == 0 ? '+' : '-';
        start_t = clock();
        SingleEndMapping(read_seqs[j], genome, hash_table, strand, AG_WILDCARD,
                         seed_len, map_results[j]);
        sum_t += clock() - start_t;
      }
      fprintf(stderr, "[%.3lf SECONDS MAPPING TIME PASSED]\n",
              static_cast<double>(sum_t / CLOCKS_PER_SEC));
    }
    OutputSingleEndResults(fout, map_results, num_of_reads, read_names,
                           read_seqs, read_scores, genome);

    if (num_of_reads < n_reads_to_process)
      break;
  }
  fclose(fin);
  fclose(fout);

  fprintf(stderr, "[MAPPING TAKES %.3lf SECONDS]\n",
          static_cast<double>(sum_t / CLOCKS_PER_SEC));
}
