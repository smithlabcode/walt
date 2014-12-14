/*
 * This is the main function for building index for reference genome.
 */
#include <string>
#include <vector>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "reference.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;

int main(int argc, const char **argv) {
  try {
    string chrom_file;
    string outfile;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "build index for reference genome", "");
    opt_parse.add_opt(
        "chrom",
        'c',
        "chromosomes in FASTA file or dir \
        (the suffix of the chromosome file should be '.fa')",
        true, chrom_file);
    opt_parse.add_opt(
        "output", 'o',
        "output file name (the suffix of the file should be '.dbindex')", true,
        outfile);
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
    if (!is_valid_filename(outfile, "dbindex")) {
      cerr << "The suffix of the output file should be '.dbindex' " << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    //////////////////////////////////////////////////////////////
    // BUILD THE INDEX
    //
    vector<string> chrom_files;
    IdentifyChromosomes(chrom_file, chrom_files);

    Genome genome;
    HashTable hash_table;
    TIME_INFO(ReadGenome(chrom_files, &genome), "READ CHROMOSOMES");
    TIME_INFO(BuildHashTable(genome, &hash_table), "BUILD HASH TABLE");
    TIME_INFO(WriteIndex(outfile, genome, hash_table), "WRITE INDEX");
    //TIME_INFO(TestHashTable(genome, hash_table), "TEST HASH TABLE");
  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
