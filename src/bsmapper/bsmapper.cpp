/*
 * This is the main function for bsmapper.
 */

#include <omp.h>

#include <fstream>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "mapping.hpp"
#include "reference.hpp"

using std::cerr;
using std::endl;
using std::ofstream;

/* load reads from reads file, each time load n_reads_to_process reads,
 * start from  read_start_idx */
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
      uint32_t space_pos = line.find_first_of(' ');
      if (space_pos == string::npos) {
        read_names.push_back(line.substr(1));
      } else {
        read_names.push_back(line.substr(1, space_pos - 1));
      }
    }
    ++line_count;
  }
}

int main(int argc, const char **argv) {
  try {
    string reads_file;
    string index_file;
    string outfile;
    size_t max_mismatches = std::numeric_limits<size_t>::max();
    size_t n_reads_to_process = std::numeric_limits<size_t>::max();
    int num_top_diags = 50;
    int num_of_threads = 1;

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
        "reads",
        'r',
        "reads file (the suffix of the reads file should be '.fastq', '.fq', '.fasta' or '.fa')",
        true, reads_file);

    opt_parse.add_opt("output", 'o', "output file name", true, outfile);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", false,
                      max_mismatches);
    opt_parse.add_opt("number", 'N', "number of reads to map at one loop",
                      false, n_reads_to_process);
    opt_parse.add_opt("diag", 'd', "number of diags for linear check", false,
                      num_top_diags);
    opt_parse.add_opt("thread", 't', "number of threads", false,
                      num_of_threads);
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
      cerr
          << "The suffix of the reads file should be '.fastq', '.fq', '.fasta' or '.fa'"
          << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    //////////////////////////////////////////////////////////////
    // LOAD THE INDEX
    Genome genome;
    TIME_INFO(ReadIndex(index_file, &genome), "READ INDEX");

    //////////////////////////////////////////////////////////////
    // SET NUMBER OF THREADS
    int num_of_threads_in_system = omp_get_max_threads();
    cerr << "[THERE ARE TOTAL " << num_of_threads_in_system
         << " THREADS IN SYSMTE]" << endl;
    if (num_of_threads_in_system < num_of_threads) {
      cerr << "[THERE ARE " << num_of_threads_in_system << " THREADS IN SYSMTE]"
           << endl;
      num_of_threads = num_of_threads_in_system;
    }
    omp_set_dynamic(1);
    omp_set_num_threads(num_of_threads);
    cerr << "[" << num_of_threads << " THREADS FOR MAPPING]" << endl;
    //////////////////////////////////////////////////////////////
    // LOAD THE READS
    vector<string> read_names;
    vector<string> read_seqs;
    clock_t start_t, end_t;
    start_t = clock();
    ofstream fout(outfile.c_str());
    for (uint64_t i = 0;; i += n_reads_to_process) {
      LoadReadsFromFastqFile(reads_file, i, n_reads_to_process, read_names,
                             read_seqs);
      uint32_t num_of_reads = read_seqs.size();
      if (num_of_reads == 0)
        break;
      uint32_t read_width = read_seqs[0].size();
      if (max_mismatches == std::numeric_limits<size_t>::max()) {
        max_mismatches = static_cast<size_t>(0.07 * read_width);
      }

      BestMatch best_match(0, 0, 0, max_mismatches);
      vector<BestMatch> map_results(num_of_reads, best_match);
      cerr << "[START MAPPING]" << endl;
#pragma omp parallel for
      for (uint32_t j = 0; j < num_of_reads; ++j) {
        SingleEndMapping(read_seqs[j].c_str(), genome, map_results[j]);
      }
#pragma omp barrier

      for (uint32_t j = 0; j < num_of_reads; ++j) {
        fout << genome[map_results[j].chrom_id].name << "\t"
             << map_results[j].chrom_pos << "\t" << read_names[j] << "\t"
             << map_results[j].mismatch << "\t"
             << genome[map_results[j].chrom_id].strand << endl;
      }
      if (read_seqs.size() < n_reads_to_process)
        break;
    }
    fout.close();

    end_t = clock();
    printf("[MAPPING TAKES %.3lf SECONDS]\n",
           (double) ((end_t - start_t) / CLOCKS_PER_SEC));

  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
