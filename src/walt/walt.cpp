/*
 *    This is the main function for bisulfite sequencing mapping.
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

#include <vector>
#include <string>
#include <omp.h>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

#include "paired.hpp"
#include "mapping.hpp"
#include "reference.hpp"

int main(int argc, const char **argv) {
  srand (time(NULL));
  try {
    string command = argv[0];
    bool help_info = false;
    for (int i = 1; i < argc; i++) {
      command += " ";
      command += argv[i];
      if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "-about") == 0
          || strcmp(argv[i], "-?") == 0) {
        help_info = true;
      }
    }

    if (argc > 1 && help_info == false) {
      /* show the command line one the screen */
      fprintf(stderr, "[WELCOME TO WALT v%s]\n", walt_version);
      fprintf(stderr, "[%s", argv[0]);
      for (int i = 1; i < argc; i++) {
        fprintf(stderr, " %s", argv[i]);
      }
      fprintf(stderr, "]\n");
    }

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

    /* WALT supports Tab-delimited SAM and MR output formats.
     * By default, WALT produces SAM format output files. To
     * get MR format, the suffix of the output file should be
     * ".mr". */
    bool SAM = true;

    /* trimming adaptor sequence */
    string adaptor;

    /* paired-end or single-end mapping */
    bool has_paired_end_reads = false;

    /* output ambiguous or unmapped reads or not, both are false by default */
    bool ambiguous = false;
    bool unmapped = false;

    /* AG_WILDCARD is false by default, which means that all Cs
     * in the reads and genome are converted to Ts.
     * If AG_WILDCARD is true, all Gs in the reads and genome
     * are converted to As. If option is only for single-end mapping */
    bool AG_WILDCARD = false;

    /* post-bisulfite adaptor tagging  */
    bool PBAT = false;

    /* maximum allowed mismatches */
    uint32_t max_mismatches = 6;

    /* number of reads to map at one loop */
    uint32_t n_reads_to_process = 5000000;
    const uint32_t PROCESS_LIMIT = 10000000;

    /* ignore seeds which have more than b candidates */
    uint32_t b = 5000;

    /* paired-end reads: keep top k genome positions for each in the pair */
    uint32_t top_k = 50;

    /* max fragment length for paired end reads */
    int frag_range = 1000;

    /* number of threads for mapping */
    int num_of_threads = 1;
    FILE * fambiguous = NULL, * funmapped = NULL, * fout, * fstat;

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
    opt_parse.add_opt("ag-wild", 'A',
                      "map using A/G bisulfite wildcards (single-end)", false,
                      AG_WILDCARD);
    opt_parse.add_opt("pbat", 'P',
                      "map post-bisulfite adaptor tagging reads", false,
                      PBAT);
    opt_parse.add_opt("bucket", 'b', "maximum candidates for a seed",
                      false, top_k);
    opt_parse.add_opt("topk", 'k',
                      "maximum allowed mappings for a read (paired-end)", false,
                      top_k);
    opt_parse.add_opt("fraglen", 'L', "max fragment length (paired-end)", false,
                      frag_range);
    opt_parse.add_opt("thread", 't', "number of threads for mapping", false,
                      num_of_threads);

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

    if (!leftover_args.empty()) {
      fprintf(stderr, "%s\n", opt_parse.help_message().c_str());
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    bool get_empty_fields = false;
    if (!reads_file_s.empty()) {
      v_reads_file_s = smithlab::split(reads_file_s, ",", get_empty_fields);
      for (uint32_t i = 0; i < v_reads_file_s.size(); ++i) {
        if (!is_valid_filename(v_reads_file_s[i], "fastq")
            && !is_valid_filename(v_reads_file_s[i], "fq")) {
          fprintf(stderr,
                  "The suffix of the reads file should be '.fastq', '.fq'\n");
          return EXIT_FAILURE;
        }
      }
    } 
    if (!reads_file_p1.empty() && !reads_file_p2.empty()) {
      has_paired_end_reads = true;
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
    else if (reads_file_p1.empty() ^ reads_file_p2.empty()){
      fprintf(stderr, "When using -1 and -2 for pair-end reads, both of them \n\
          should be specified.\n");
      return EXIT_FAILURE;
    }

    // CHECK OPTIONS
    omp_set_dynamic(0);
    omp_set_num_threads(num_of_threads);
    fprintf(stderr, "[MAXIMUM NUMBER OF MISMATCHES IS %u]\n", max_mismatches);
    fprintf(stderr, "[NUMBER OF THREADS FOR MAPPING IS %d]\n", num_of_threads);

    if (n_reads_to_process > PROCESS_LIMIT) {
      n_reads_to_process = PROCESS_LIMIT;
    }

    if (has_paired_end_reads && top_k < 2) {
      fprintf(stderr, "-k option should be at least 2 for paired-end reads\n");
      return EXIT_FAILURE;
    }

    if (has_paired_end_reads && top_k > 300) {
      fprintf(stderr,
              "-k option should be less than 300 for paired-end reads\n");
      return EXIT_FAILURE;
    }

    size_t suffix_pos = output_file.find_last_of(".");
    if (suffix_pos == string::npos) {
      SAM = true;
    } else if (".mr" == output_file.substr(suffix_pos)) {
      SAM = false;
    }
    fout = fopen(output_file.c_str(), "w");
    if (!fout) {
      throw SMITHLABException("cannot open input file " + output_file);
    }
    fprintf(stderr, "[OUTPUT MAPPING RESULTS TO %s]\n", output_file.c_str());
    if(SAM) {
      SAMHead(index_file, command, fout);
    }
    if (ambiguous && !SAM) {
      fambiguous = fopen(string(output_file + ".ambiguous").c_str(), "w");
      if (!fambiguous) {
        throw SMITHLABException(
            "cannot open output file " + string(output_file + ".ambiguous"));
      }
    }
    if (unmapped && !SAM) {
      funmapped = fopen(string(output_file + ".unmapped").c_str(), "w");
      if (!funmapped) {
        throw SMITHLABException(
            "cannot open output file " + string(output_file + ".unmapped"));
      }
    }
    fstat = fopen(string(output_file.substr(0, suffix_pos) \
          + ".mapstats").c_str(), "w");

    ShowGenomeInfo(index_file);
    // To-do: Load the mapping index here. There is no need to load index repeatedly.
    //
    // Mapping
    if (!reads_file_s.empty()) {
      for (uint32_t i = 0; i < v_reads_file_s.size(); ++i) {
        ProcessSingledEndReads(index_file, v_reads_file_s[i],
                               fout, fstat, n_reads_to_process,
                               max_mismatches, b, adaptor, PBAT, AG_WILDCARD,
                               fambiguous, funmapped, SAM);
      }
    }
    if (has_paired_end_reads){
      for (uint32_t i = 0; i < v_reads_file_p1.size(); ++i) {
        ProcessPairedEndReads(index_file, v_reads_file_p1[i],
                              v_reads_file_p2[i], fout, fstat,
                              n_reads_to_process, max_mismatches, b, adaptor,
                              PBAT, top_k, frag_range, fambiguous, funmapped,
                              SAM);
      }
    }
    fclose(fout);
    fclose(fstat);
    if (fambiguous) fclose(fambiguous);
    if (funmapped) fclose(funmapped);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
