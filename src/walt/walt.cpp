/*
 * This is the main function for bisulfite sequencing mapping.
 * Copyright [2015] <Andrew D. Smith, Ting Chen>
 */

#include <vector>
#include <string>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

#include "paired.hpp"
#include "mapping.hpp"
#include "reference.hpp"

int main(int argc, const char **argv) {
  try {
    fprintf(stderr, "[WELCOME TO WALT v0.0]\n");
    fprintf(stderr, "[%s", argv[0]);
    for (int i = 1; i < argc; i++) {
      fprintf(stderr, " %s", argv[i]);
    }
    fprintf(stderr, "]\n");

    /* singled-end reads file, comma-separated list of files */
    string reads_file_s;
    vector<string> v_reads_file_s;

    /* paired-end reads files, comma-separated list of files*/
    string reads_file_p1;
    string reads_file_p2;
    vector<string> v_reads_file_p1;
    vector<string> v_reads_file_p2;

    /* index file*/
    string index_file;

    /* output file */
    string output_file;
    vector<string> v_output_file;

    /* output SAM format */
    bool SAM = false;

    /* trimming adaptor sequence */
    string adaptor;

    bool is_paired_end_reads = false;
    bool AG_WILDCARD = false;
    bool ambiguous = false;
    bool unmapped = false;

    uint32_t max_mismatches = MAX_UINT32;
    uint32_t n_reads_to_process = MAX_UINT32;

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
        "reads",
        'r',
        "comma-separated list of read files for singled-end mapping \
         (the suffix of read files should be '.fastq' or '.fq')",
        false, reads_file_s);
    opt_parse.add_opt(
        "reads1",
        '1',
        "comma-separated list of read files for mate 1 \
         (the suffix of read files should be '.fastq' or '.fq')",
        false, reads_file_p1);
    opt_parse.add_opt(
        "reads2",
        '2',
        "comma-separated list of read files for mate 2 \
        (the suffix of read files should be '.fastq' or '.fq')",
        false, reads_file_p2);
    opt_parse.add_opt("output", 'o', "output file name", true, output_file);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", false,
                      max_mismatches);
    opt_parse.add_opt("number", 'N', "number of reads to map at one loop",
                      false, n_reads_to_process);
    opt_parse.add_opt(
        "ambiguous",
        'a',
        "randomly output one mapped position for ambiguous \
        reads in a separated file",
        false, ambiguous);
    opt_parse.add_opt("unmapped", 'u',
                      "output unmapped reads in a separated file", false,
                      unmapped);
    opt_parse.add_opt("clip", 'C', "clip the specified adaptor", false,
                      adaptor);
    opt_parse.add_opt("ag-wild", 'A', "map using A/G bisulfite wildcards",
                      false, AG_WILDCARD);
    opt_parse.add_opt("topk", 'k',
                      "maximum allowed mappings for a read (paired-end)", false,
                      top_k);
    opt_parse.add_opt("fraglen", 'L', "max fragment length (paired-end)", false,
                      frag_range);

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

    bool get_empty_fields = false;
    if (!is_paired_end_reads) {
      v_reads_file_s = smithlab::split(reads_file_s, ",", get_empty_fields);
      for (uint32_t i = 0; i < v_reads_file_s.size(); ++i) {
        if (!is_valid_filename(v_reads_file_s[i], "fastq")
            && !is_valid_filename(v_reads_file_s[i], "fq")) {
          fprintf(stderr,
                  "The suffix of the reads file should be '.fastq', '.fq'\n");
          return EXIT_FAILURE;
        }
      }
    } else {
      v_reads_file_p1 = smithlab::split(reads_file_p1, ",", get_empty_fields);
      v_reads_file_p2 = smithlab::split(reads_file_p2, ",", get_empty_fields);
      if (v_reads_file_p1.size() != v_reads_file_p2.size()) {
        fprintf(
            stderr,
            "For paired-end mapping, mate 1 and mate 2 should \n\
                have the same number of files, and the paired files \n\
                should be in the same order.");
        return EXIT_FAILURE;
      }
      for (uint32_t i = 0; i < v_reads_file_p1.size(); ++i) {
        if (!is_valid_filename(v_reads_file_p1[i], "fastq")
            && !is_valid_filename(v_reads_file_p1[i], "fq")) {
          fprintf(stderr,
                  "The suffix of the reads file should be '.fastq', '.fq'\n");
          return EXIT_FAILURE;
        }
      }
    }

    if (!is_paired_end_reads) {
      if (v_reads_file_s.size() == 1) {
        v_output_file.push_back(output_file);
      } else {
        char output_filename[1000];
        for (uint32_t i = 0; i < v_reads_file_s.size(); ++i) {
          sprintf(output_filename, "%s_s%u", output_file.c_str(), i);
          v_output_file.push_back(output_filename);
        }
      }
    } else {
      if (v_reads_file_p1.size() == 1) {
        v_output_file.push_back(output_file);
      } else {
        char output_filename[1000];
        for (uint32_t i = 0; i < v_reads_file_p1.size(); ++i) {
          sprintf(output_filename, "%s_p%u", output_file.c_str(), i);
          v_output_file.push_back(output_filename);
        }
      }
    }

    for (uint32_t i = 0; i < v_reads_file_p2.size(); ++i) {
      if (!is_valid_filename(v_reads_file_p2[i], "fastq")
          && !is_valid_filename(v_reads_file_p2[i], "fq")) {
        fprintf(stderr,
                "The suffix of the reads file should be '.fastq', '.fq'\n");
        return EXIT_FAILURE;
      }
    }

    uint32_t suffix_pos = output_file.find_last_of(".");
    if (".sam" == output_file.substr(suffix_pos)) {
      SAM = true;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    //////////////////////////////////////////////////////////////
    // CHECK OPTIONS
    if (max_mismatches == MAX_UINT32) {
      max_mismatches = 6;
      fprintf(stderr, "[MAXIMUM NUMBER OF MISMATCHES IS %u]\n", max_mismatches);
    }

    if (n_reads_to_process > 5000000) {
      n_reads_to_process = 5000000;
    }

    if (is_paired_end_reads && top_k < 2) {
      fprintf(stderr, "-k option should be at least 2 for paired-end reads\n");
      return EXIT_FAILURE;
    }

    if (is_paired_end_reads && top_k > 300) {
      fprintf(stderr,
              "-k option should be less than 300 for paired-end reads\n");
      return EXIT_FAILURE;
    }
    ShowGenomeInfo(index_file);

    //////////////////////////////////////////////////////////////
    // Mapping
    if (!is_paired_end_reads) {
      for (uint32_t i = 0; i < v_reads_file_s.size(); ++i) {
        ProcessSingledEndReads(index_file, v_reads_file_s[i], v_output_file[i],
                               n_reads_to_process, max_mismatches, adaptor,
                               AG_WILDCARD, ambiguous, unmapped, SAM);
      }
    } else {
      for (uint32_t i = 0; i < v_reads_file_p1.size(); ++i) {
        ProcessPairedEndReads(index_file, v_reads_file_p1[i],
                              v_reads_file_p2[i], v_output_file[i],
                              n_reads_to_process, max_mismatches, adaptor,
                              top_k, frag_range, ambiguous, unmapped, SAM);
      }
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
