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

  cerr << "[LOADING READS]" << endl;

  string line;
  int line_code = 0;
  uint64_t line_count = 0;
  num_of_reads = 0;
  uint64_t lim = n_reads_to_process * 4;
  while (line_count < lim && getline(fin, line)) {
    switch (line_code) {
      case 0: {
        uint32_t space_pos = line.find_first_of(' ');
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

int main(int argc, const char **argv) {
  try {
    string reads_file;
    string index_file;
    string outfile;
    size_t max_mismatches = std::numeric_limits<size_t>::max();
    size_t n_reads_to_process = std::numeric_limits<size_t>::max();

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
        true, reads_file);

    opt_parse.add_opt("output", 'o', "output file name", true, outfile);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", false,
                      max_mismatches);
    opt_parse.add_opt("number", 'N', "number of reads to map at one loop",
                      false, n_reads_to_process);

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
      cerr << "The suffix of the index file should be '.dbindex' " << endl;
      return EXIT_SUCCESS;
    }
    if (!is_valid_filename(reads_file, "fastq")
        && !is_valid_filename(reads_file, "fq")) {
      cerr << "The suffix of the reads file should be '.fastq', '.fq'" << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    //////////////////////////////////////////////////////////////
    // LOAD THE INDEX
    Genome genome;
    HashTable hash_table;
    TIME_INFO(ReadIndex(index_file, &genome, &hash_table), "READ INDEX");

    //////////////////////////////////////////////////////////////
    // LOAD THE READS
    if (n_reads_to_process > 1000000) {
      n_reads_to_process = 1000000;
    }
    vector<string> read_names(n_reads_to_process);
    vector<string> read_seqs(n_reads_to_process);
    vector<string> read_scores(n_reads_to_process);

    clock_t start_t, end_t;
    start_t = clock();

    ifstream fin(reads_file.c_str());
    if (!fin) {
      throw SMITHLABException("cannot open input file " + reads_file);
    }
    ofstream fout(outfile.c_str());
    uint32_t num_of_reads;
    vector<BestMatch> map_results(n_reads_to_process);
    for (uint64_t i = 0;; i += n_reads_to_process) {
      LoadReadsFromFastqFile(fin, i, n_reads_to_process, num_of_reads,
                             read_names, read_seqs, read_scores);
      if (num_of_reads == 0)
        break;
      uint32_t read_width = read_seqs[0].size();
      if (max_mismatches == std::numeric_limits<size_t>::max()) {
        max_mismatches = static_cast<size_t>(0.07 * read_width);
        cerr << "[MAXIMUM NUMBER OF MISMATCHES IS " << max_mismatches << "]"
             << endl;
      }

      BestMatch best_match(0, 0, 0, max_mismatches);
      for (uint32_t j = 0; j < num_of_reads; ++j) {
        map_results[j] = best_match;
      }

      cerr << "[START MAPPING READS FROM " << i << " TO " << num_of_reads + i
           << "]" << endl;

      MappingAllReads(read_seqs, num_of_reads, genome, hash_table, map_results);

      for (uint32_t j = 0; j < num_of_reads; ++j) {
        if (map_results[j].times == 0 || map_results[j].times > 1)
          continue;
        uint32_t start_pos = map_results[j].chrom_pos;
        if ('-' == genome[map_results[j].chrom_id].strand) {
          start_pos = genome[map_results[j].chrom_id].length
              - map_results[j].chrom_pos - read_seqs[j].size();
        }
        uint32_t end_pos = start_pos + read_seqs[j].size();

        fout << genome[map_results[j].chrom_id].name << "\t" << start_pos
             << "\t" << end_pos << "\t" << read_names[j] << "\t"
             << map_results[j].mismatch << "\t"
             << genome[map_results[j].chrom_id].strand << "\t" << read_seqs[j]
             << "\t" << read_scores[j] << endl;
      }

      if (num_of_reads < n_reads_to_process)
        break;
    }
    fin.close();
    fout.close();

    end_t = clock();
    fprintf(stderr, "[MAPPING TAKES %.3lf SECONDS]\n",
            static_cast<double>((end_t - start_t) / CLOCKS_PER_SEC));
  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
