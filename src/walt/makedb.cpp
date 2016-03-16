/*
 *    This is the main function to build index for reference genome.
 *
 *    Copyright (C) 2015 University of Southern California
 *                       Andrew D. Smith and Ting Chen
 *
 *    Authors: Haifeng Chen, Andrew D. Smith and Ting Chen
 *
 *    This file is part of WALT.
 *
 *    WALT is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    WALT is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with WALT.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <set>
#include <string>
#include <vector>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "reference.hpp"

void BuildIndex(const vector<string>& chrom_files, const int& indicator,
                const string& output_file, uint32_t& size_of_index) {
  switch (indicator) {
    case 0:
      fprintf(stderr, "[BIULD INDEX FOR FORWARD STRAND (C->T)]\n");
      break;
    case 1:
      fprintf(stderr, "[BIULD INDEX FOR REVERSE STRAND (C->T)]\n");
      break;
    case 2:
      fprintf(stderr, "[BIULD INDEX FOR FORWARD STRAND (G->A)]\n");
      break;
    case 3:
      fprintf(stderr, "[BIULD INDEX FOR REVERSE STRAND (G->A)]\n");
  }

  Genome genome;
  HashTable hash_table;
  ReadGenome(chrom_files, genome);

  if (indicator % 2) {
    ReverseComplementGenome(genome);
  }

  if (indicator == 0 || indicator == 1) {
    C2T(genome.sequence);
  } else {
    G2A(genome.sequence);
  }

  set<uint32_t> extremal_large_bucket;
  CountBucketSize(genome, hash_table, extremal_large_bucket);
  HashToBucket(genome, hash_table, extremal_large_bucket);
  SortHashTableBucket(genome, hash_table);
  WriteIndex(output_file, genome, hash_table);

  size_of_index =
      hash_table.index_size > size_of_index ?
          hash_table.index_size : size_of_index;
}

int main(int argc, const char **argv) {
  srand (time(NULL));
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
      fprintf(stderr, "%s\n", opt_parse.help_message().c_str());
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      fprintf(stderr, "%s\n", opt_parse.about_message().c_str());
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      fprintf(stderr, "%s\n", opt_parse.option_missing_message().c_str());
      return EXIT_SUCCESS;
    }
    if (!is_valid_filename(outfile, "dbindex")) {
      fprintf(stderr, "The suffix of the output file should be '.dbindex'\n");
      return EXIT_FAILURE;
    }
    if (outfile.size() > 1000) {
      fprintf(stderr, "The output file name is too long, "
              "please select a shorter name\n");
      return EXIT_FAILURE;;
    }
    if (!leftover_args.empty()) {
      fprintf(stderr, "%s\n", opt_parse.help_message().c_str());
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
    ////////// BUILD INDEX FOR FORWARD STRAND (C->T)
    BuildIndex(chrom_files, 0, outfile + "_CT00", size_of_index);

    ////////// BUILD INDEX FOR REVERSE STRAND (C->T)
    BuildIndex(chrom_files, 1, outfile + "_CT01", size_of_index);

    ////////// BUILD INDEX FOR FORWARD STRAND (G->A)
    BuildIndex(chrom_files, 2, outfile + "_GA10", size_of_index);

    ////////// BUILD INDEX FOR REVERSE STRAND (G->A)
    BuildIndex(chrom_files, 3, outfile + "_GA11", size_of_index);

    Genome genome;
    ReadGenome(chrom_files, genome);
    WriteIndexHeadInfo(outfile, genome, size_of_index);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
