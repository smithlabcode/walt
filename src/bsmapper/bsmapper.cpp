/*
 * This is the main function for bisulfite sequence mapping.
 * Copyright [2015] < >
 */

#include <vector>
#include <string>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "paired.hpp"
#include "mapping.hpp"
#include "reference.hpp"

/* get the length of the reads in the fastq file */
uint32_t GetReadLength(const string& reads_file) {
  FILE * fin = fopen(reads_file.c_str(), "r");
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
      fprintf(stderr, "All the reads should have the same length.\n"
              "Please check the reads file.\n");
      return EXIT_FAILURE;
    }
  }

  fclose(fin);

  return read_len;
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

    uint32_t max_mismatches = MAX_UINT32;
    uint32_t n_reads_to_process = MAX_UINT32;
    uint32_t seed_len = MAX_UINT32;

    /* paired-end reads: keep top k genome positions for each in the pair */
    uint32_t top_k = 50;

    /* max fragment length for paired end reads */
    int frag_range = 1000;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "map Illumina BS-seq reads",
                           "");
    opt_parse.add_opt(
        "index",
        'i',
        "index file created by makedb command \
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
        "reads2", '2',
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
                      top_k);
    opt_parse.add_opt("fraglen", 'L', "max fragment length", false, frag_range);

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

    if (!is_valid_filename(index_file, "dbindex")) {
      fprintf(stderr, "The suffix of the index file should be '.dbindex'\n");
      return EXIT_FAILURE;
    }

    if (!reads_file_s.empty() && reads_file_p1.empty()
        && reads_file_p2.empty()) {
      is_paired_end_reads = false;
    } else if (reads_file_s.empty() && !reads_file_p1.empty()
        && !reads_file_p2.empty()) {
      is_paired_end_reads = true;
    } else {
      fprintf(stderr, "Please use -r option to set singled-end reads, \n"
              "-1 and -2 options to set paired-end reads\n");
      return EXIT_FAILURE;
    }

    if (!is_paired_end_reads && !is_valid_filename(reads_file_s, "fastq")
        && !is_valid_filename(reads_file_s, "fq")) {
      fprintf(stderr,
              "The suffix of the reads file should be '.fastq', '.fq'\n");
      return EXIT_FAILURE;
    }
    if (is_paired_end_reads) {
      if ((!is_valid_filename(reads_file_p1, "fastq")
          && !is_valid_filename(reads_file_p1, "fq"))
          || (!is_valid_filename(reads_file_p1, "fastq")
              && !is_valid_filename(reads_file_p1, "fq"))) {
        fprintf(stderr,
                "The suffix of the reads file should be '.fastq', '.fq'\n");
        return EXIT_FAILURE;
      }
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    //////////////////////////////////////////////////////////////
    // CHECK OPTIONS
    uint32_t read_len;
    if (!is_paired_end_reads) {
      read_len = GetReadLength(reads_file_s);
    } else {
      uint32_t read_len1 = GetReadLength(reads_file_p1);
      uint32_t read_len2 = GetReadLength(reads_file_p2);
      if (read_len1 != read_len2) {
        fprintf(stderr, "All the reads should have the same length. \n"
                "Please check the reads file.\n");
        return EXIT_FAILURE;
      }
      read_len = read_len1;
    }
    fprintf(stderr, "[READ LENGTH IS %u]\n", read_len);

    if (read_len < HASHLEN) {
      fprintf(stderr, "The length of the reads should be at least %u\n",
              HASHLEN);
      return EXIT_FAILURE;
    }

    if (max_mismatches == MAX_UINT32) {
      max_mismatches = static_cast<size_t>(0.07 * read_len);
      fprintf(stderr, "[MAXIMUM NUMBER OF MISMATCHES IS %u]\n", max_mismatches);
    }

    if (seed_len == MAX_UINT32) {
      seed_len = (read_len - SEEPATTERNLEN + 1) / SEEPATTERNLEN;
    }

    if (seed_len < F2SEEDWIGTH) {
      fprintf(stderr, "The seed length should be at least %u\n", F2SEEDWIGTH);
      return EXIT_FAILURE;
    }

    if (seed_len > F2SEEDPOSITION_SIZE) {
      fprintf(stderr, "The seed length should be no more than %u\n",
              F2SEEDPOSITION_SIZE);
      return EXIT_FAILURE;
    }

    if (F2SEEDPOSITION[seed_len - 1] >= read_len - SEEPATTERNLEN) {
      fprintf(stderr,
              "[THE SEED LENGTH SHOULD BE SHORTER FOR THIS READ LENGTH]\n");
      return EXIT_FAILURE;
    } else {
      fprintf(stderr, "[SEED LENGTH IS %u]\n", seed_len);
    }

    if (n_reads_to_process > 5000000) {
      n_reads_to_process = 5000000;
    }

    if (is_paired_end_reads && top_k < 2) {
      fprintf(stderr, "-k option should be at least 2 for paired-end reads\n");
      return EXIT_FAILURE;
    }

    if (is_paired_end_reads && top_k > 300) {
      fprintf(stderr, "-k option should be less than 300 for paired-end reads\n");
      return EXIT_FAILURE;
    }

    //////////////////////////////////////////////////////////////
    // Mapping
    if (!is_paired_end_reads) {
      ProcessSingledEndReads(index_file, reads_file_s, output_file,
                             n_reads_to_process, max_mismatches, read_len,
                             seed_len, AG_WILDCARD);
    } else {
      ProcessPairedEndReads(index_file, reads_file_p1, reads_file_p2,
                            output_file, n_reads_to_process, max_mismatches,
                            read_len, seed_len, top_k, frag_range);
    }
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
