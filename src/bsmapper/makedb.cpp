/*
 * This is the main function for building index for reference genome.
 * Copyright [2014] < >
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

void BuildIndex(const vector<string>& chrom_files, const int& indicator,
                const string& output_file, uint32_t& size_of_index) {
  switch (indicator) {
    case 0:
      cerr << "[BIULD INDEX FOR FORWARD STRAND (C->T)]" << endl;
      break;
    case 1:
      cerr << "[BIULD INDEX FOR REVERSE STRAND (C->T)]" << endl;
      break;
    case 2:
      cerr << "[BIULD INDEX FOR FORWARD STRAND (A->G)]" << endl;
      break;
    case 3:
      cerr << "[BIULD INDEX FOR REVERSE STRAND (A->G)]" << endl;
  }

  Genome genome;
  HashTable hash_table;
  ReadGenome(chrom_files, genome);

  if (indicator % 2) {
    ReverseGenome(genome);
  }

  if (indicator == 0 || indicator == 1) {
    C2T(genome.sequence);
  } else {
    A2G(genome.sequence);
  }

  CountBucketSize(genome, hash_table);
  HashToBucket(genome, hash_table);
  SortHashTableBucket(genome, hash_table);
  WriteIndex(output_file, genome, hash_table);

  size_of_index =
      hash_table.index_size > size_of_index ?
          hash_table.index_size : size_of_index;
}

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
    if (outfile.size() > 1000) {
      cerr << "The output file name is too long, please select a shorter name"
          << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    //////////////////////////////////////////////////////////////
    // READ GENOME
    //
    vector<string> chrom_files;
    IdentifyChromosomes(chrom_file, chrom_files);

    //////////////////////////////////////////////////////////////
    // BUILD  INDEX
    //
    uint32_t size_of_index = 0;
    //////////BUILD INDEX FOR FORWARD STRAND (C->T)
    BuildIndex(chrom_files, 0, outfile + "_CT00", size_of_index);

    //////////BUILD INDEX FOR REVERSE STRAND (C->T)
    BuildIndex(chrom_files, 1, outfile + "_CT01", size_of_index);

    //////////BUILD INDEX FOR FORWARD STRAND (A->G)
    BuildIndex(chrom_files, 2, outfile + "_AG10", size_of_index);

    //////////BUILD INDEX FOR REVERSE STRAND (A->G)
    BuildIndex(chrom_files, 3, outfile + "_AG11", size_of_index);

    Genome genome;
    ReadGenome(chrom_files, genome);
    WriteIndexHeadInfo(outfile, genome, size_of_index);

  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
