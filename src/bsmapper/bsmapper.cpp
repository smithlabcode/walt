/*
 * This is the main function for bsmapper.
 * Copyright [2014] < >
 */

#include <vector>
#include <string>
#include <limits>
#include <fstream>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "mapping.hpp"
#include "reference.hpp"

using std::cerr;
using std::endl;
using std::vector;
using std::ofstream;
using std::ifstream;

/* load reads from reads file, each time load n_reads_to_process reads,
 * start from  read_start_idx */
void LoadReadsFromFastqFile(ifstream &fin, const uint64_t read_start_idx,
                            const uint64_t n_reads_to_process,
                            uint32_t& num_of_reads, vector<string>& read_names,
                            vector<string>& read_seqs,
                            vector<string>& read_scores) {
  string line;
  int line_code = 0;
  uint64_t line_count = 0;
  num_of_reads = 0;
  uint64_t lim = n_reads_to_process * 4;
  while (line_count < lim && getline(fin, line)) {
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

/* get the length of the reads in the fastq file */
uint32_t GetReadLength(const string& reads_file) {
  ifstream fin(reads_file.c_str());
  if (!fin) {
    throw SMITHLABException("cannot open input file " + reads_file);
  }
  srand(time(NULL));
  int rand_num_reads = rand() % 100 + 10;
  vector<string> read_names(rand_num_reads);
  vector<string> read_seqs(rand_num_reads);
  vector<string> read_scores(rand_num_reads);
  uint32_t num_of_reads = 0;
  LoadReadsFromFastqFile(fin, 0, rand_num_reads, num_of_reads, read_names,
                         read_seqs, read_scores);
  if (read_seqs.size() == 0)
    return 0;
  uint32_t read_len = read_seqs[0].size();
  for (uint32_t i = 1; i < num_of_reads; ++i) {
    if (read_seqs[i].size() != read_len) {
      cerr << "All the reads should have the same length. "
          << "Please check the reads file." << endl;
      return EXIT_FAILURE;
    }
  }

  fin.close();

  return read_len;
}

void OutputSingleEndResults(ofstream& fout,
                            const vector<BestMatch>& map_results,
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

    fout << genome.name[chr_id] << "\t" << start_pos << "\t" << end_pos << "\t"
        << read_names[j] << "\t" << map_results[j].mismatch << "\t"
        << map_results[j].strand << "\t" << read_seqs[j] << "\t"
        << read_scores[j] << endl;
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
  ifstream fin(reads_file_s.c_str());
  if (!fin) {
    throw SMITHLABException("cannot open input file " + reads_file_s);
  }
  clock_t start_t;
  uint64_t sum_t = 0;

  ofstream fout(output_file.c_str());
  uint32_t num_of_reads;
  for (uint64_t i = 0;; i += n_reads_to_process) {
    LoadReadsFromFastqFile(fin, i, n_reads_to_process, num_of_reads, read_names,
                           read_seqs, read_scores);
    if (num_of_reads == 0)
      break;

    //Initialize the results
    BestMatch best_match(0, 0, '+', max_mismatches);
    for (uint32_t j = 0; j < num_of_reads; ++j) {
      map_results[j] = best_match;
    }

    cerr << "[START MAPPING READS FROM " << i << " TO " << num_of_reads + i
        << "]" << endl;
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
  fin.close();
  fout.close();

  fprintf(stderr, "[MAPPING TAKES %.3lf SECONDS]\n",
          static_cast<double>(sum_t / CLOCKS_PER_SEC));
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

  vector<vector<string> > index_names(2, vector<string>(2));
  index_names[0][0] = index_file + "_CT00";
  index_names[0][1] = index_file + "_CT01";
  index_names[1][0] = index_file + "_AG10";
  index_names[1][1] = index_file + "_AG11";

  vector<vector<string> > read_names(2, vector<string>(n_reads_to_process));
  vector<vector<string> > read_seqs(2, vector<string>(n_reads_to_process));
  vector<vector<string> > read_scores(2, vector<string>(n_reads_to_process));

  vector<vector<TopCandidates> > top_results(
      2, vector<TopCandidates>(n_reads_to_process));

  vector<BestMatch> map_results(n_reads_to_process);
  vector<int> ranked_results_size(2);
  vector<vector<CandidatePosition> > ranked_results(
      2, vector<CandidatePosition>(top_k));

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

int main(int argc, const char **argv) {
  try {
    /* singled-end reads file */
    string reads_file_s;

    /* paired-end reads files*/
    string reads_file_p1;
    string reads_file_p2;

    /* index file*/
    string index_file;

    /* output file */
    string output_file;

    bool is_paired_end_reads = false;
    bool AG_WILDCARD = false;

    size_t max_mismatches = MAX_UINT32;
    size_t n_reads_to_process = MAX_UINT32;
    uint32_t seed_len = 20;

    /* paired-end reads: keep top k genome positions for each in the pair */
    uint32_t top_k = 100;
    int frag_range = 1000;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "map Illumina BS-seq reads",
                           "");
    opt_parse.add_opt(
        "index",
        'i',
        "index file created by build command \
        (the suffix of the index file should be '.dbindex')",
        true, index_file);

    opt_parse.add_opt(
        "reads", 'r',
        "reads file (the suffix of the reads file should be '.fastq' or '.fq')",
        false, reads_file_s);
    opt_parse.add_opt(
        "reads1",
        '1',
        "reads2 file (the suffix of the reads file should be '.fastq' or '.fq')",
        false, reads_file_p1);
    opt_parse.add_opt(
        "reads", '2',
        "reads file (the suffix of the reads file should be '.fastq' or '.fq')",
        false, reads_file_p2);

    opt_parse.add_opt("output", 'o', "output file name", true, output_file);
    opt_parse.add_opt("seedlength", 'l', "the length of the space seed", false,
                      seed_len);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", false,
                      max_mismatches);
    opt_parse.add_opt("number", 'N', "number of reads to map at one loop",
                      false, n_reads_to_process);

    opt_parse.add_opt("ag-wild", 'A', "map using A/G bisulfite wildcards",
                      false, AG_WILDCARD);

    opt_parse.add_opt("topk", 'k', "maximum allowed mappings for a read", false,
                      n_reads_to_process);

    opt_parse.add_opt("fraglen", 'L', "max fragment length", false, frag_range);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }

    if (!is_valid_filename(index_file, "dbindex")) {
      cerr << "The suffix of the index file should be '.dbindex'" << endl;
      return EXIT_SUCCESS;
    }

    if (!reads_file_s.empty() && reads_file_p1.empty()
        && reads_file_p2.empty()) {
      is_paired_end_reads = false;
    } else if (reads_file_s.empty() && !reads_file_p1.empty()
        && !reads_file_p2.empty()) {
      is_paired_end_reads = true;
    } else {
      cerr << "Please use -r option to set singled-end reads, "
          << "-1 and -2 options to set paired-end reads" << endl;
      return EXIT_SUCCESS;
    }

    if (!is_paired_end_reads && !is_valid_filename(reads_file_s, "fastq")
        && !is_valid_filename(reads_file_s, "fq")) {
      cerr << "The suffix of the reads file should be '.fastq', '.fq'" << endl;
      return EXIT_SUCCESS;
    }
    if (is_paired_end_reads) {
      if ((!is_valid_filename(reads_file_p1, "fastq")
          && !is_valid_filename(reads_file_p1, "fq"))
          || (!is_valid_filename(reads_file_p1, "fastq")
              && !is_valid_filename(reads_file_p1, "fq"))) {
        cerr << "The suffix of the reads file should be '.fastq', '.fq'"
            << endl;
        return EXIT_SUCCESS;
      }
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    //////////////////////////////////////////////////////////////
    // CHECK OPTIONS
    if (seed_len < F2SEEDWIGTH) {
      cerr << "The seed length should be at least " << F2SEEDWIGTH << endl;
      return EXIT_FAILURE;
    }

    if (seed_len > F2SEEDPOSITION_SIZE) {
      cerr << "The seed length should be no more than " << F2SEEDPOSITION_SIZE
          << endl;
      return EXIT_FAILURE;
    }

    uint32_t read_len;
    if (!is_paired_end_reads) {
      read_len = GetReadLength(reads_file_s);
    } else {
      uint32_t read_len1 = GetReadLength(reads_file_p1);
      uint32_t read_len2 = GetReadLength(reads_file_p2);
      if (read_len1 != read_len2) {
        cerr << "All the reads should have the same length. "
            << "Please check the reads file." << endl;
        return EXIT_FAILURE;
      }
      read_len = read_len1;
    }
    cerr << "[READ LENGTH IS " << read_len << "]" << endl;

    if (read_len < HASHLEN) {
      cerr << "The length of the reads should be at least " << HASHLEN << endl;
      return EXIT_FAILURE;
    }

    if (max_mismatches == MAX_UINT32) {
      max_mismatches = static_cast<size_t>(0.07 * read_len);
      cerr << "[MAXIMUM NUMBER OF MISMATCHES IS " << max_mismatches << "]"
          << endl;
    }

    if (F2SEEDPOSITION[seed_len - 1] >= read_len - SEEPATTERNLEN) {
      cerr << "[THE SEED LENGTH SHOULD BE SHORTER FOR THIS READ LENGTH]"
          << endl;
      return EXIT_FAILURE;
    } else {
      cerr << "[SEED LENGTH IS " << seed_len << "]" << endl;
    }

    if (n_reads_to_process > 5000000) {
      n_reads_to_process = 5000000;
    }

    if (!is_paired_end_reads) {
      top_k = 2;
    }

    if (is_paired_end_reads && top_k < 2) {
      cerr << "-k option should be at least 2 for paired-end reads" << endl;
      return EXIT_FAILURE;
    }

    if (!is_paired_end_reads) {
      ProcessSingledEndReads(index_file, n_reads_to_process, reads_file_s,
                             output_file, max_mismatches, read_len, seed_len,
                             AG_WILDCARD);
    } else {
      ProcessPairedEndReads(index_file, n_reads_to_process, reads_file_p1,
                            reads_file_p2, output_file, max_mismatches,
                            read_len, seed_len, top_k, frag_range);
    }
  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
