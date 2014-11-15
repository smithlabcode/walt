/*
 * This is the main function for bsmapper.
 */
#include <string>
#include <vector>
#include <limits>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "index.hpp"
#include "mapping.hpp"
#include "reference.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;

int main(int argc, const char **argv) {
  try {
    string reads_file;
    string index_file;
    string outfile;
    size_t max_mismatches = std::numeric_limits<size_t>::max();
    size_t n_reads_to_process = std::numeric_limits<size_t>::max();
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
    HashTable hash_table;
    ReadIndex readIndex(index_file, &genome, &hash_table);

    //////////////////////////////////////////////////////////////
    // LOAD THE READS
    vector<string> read_names;
    vector<string> read_seqs;
    for (uint64_t i = 0;; i += n_reads_to_process) {
      LoadReadsFromFastqFile(reads_file, i, n_reads_to_process, read_names,
                             read_seqs);

      uint32_t num_of_reads = read_seqs.size();
//#pragma omp parallel for
      for (uint32_t j = 0; j < num_of_reads; ++j) {

      }

      if (read_seqs.size() < n_reads_to_process)
        break;
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

