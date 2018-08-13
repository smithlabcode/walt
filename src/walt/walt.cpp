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
#include <stdexcept>

#include <sys/stat.h>

#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

#include "paired.hpp"
#include "mapping.hpp"
#include "reference.hpp"

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::to_string;

static void
split_filenames(const string &s, vector<string> &filenames) {
  std::istringstream iss(s);
  string token;
  while (std::getline(iss, token, ','))
    filenames.push_back(token);
}

static bool
valid_suffix(const string &s, const vector<string> &v) {

  bool found = false;
  for (auto && i : v)
    found = found || s.substr(s.length() - i.length()) == i;
  return found;
}


static void
validate_index_file(const string &index) {

  struct stat st;
  stat(index.c_str(), &st);
  if (!S_ISREG(st.st_mode))
    throw runtime_error("bad index file: " + index);

  vector<string> suffixes = { "_CT00",
                              "_CT01",
                              "_GA10",
                              "_GA11" };
  for (auto && i : suffixes) {
    const string suff_tab_filename = index + i;
    stat(suff_tab_filename.c_str(), &st);
    if (!S_ISREG(st.st_mode))
      throw runtime_error("bad table file: " + suff_tab_filename);
  }
}


int
main(int argc, const char **argv) {
  try {

    static const vector<string> fastq_suffixes = {".fastq", ".fq"};

    srand (time(NULL));

    string se_read_files_csv; // singled-end reads file or comma-sep list

    // paired-end reads files or comma-sep list
    string pe_read_files_end1_csv, pe_read_files_end2_csv;

    string index_file; // main index file
    string output_files_csv; // output file or comma sep list

    /* WALT supports Tab-delimited SAM and MR output formats.  By
     * default, WALT produces SAM format output files. To get MR
     * format, the suffix of the output file should be ".mr". */
    bool sam_format = false;

    string adaptor; // adaptor seq to trim at 3' end

    bool ambiguous = false; // output ambiguous reads
    bool unmapped = false; // output unmapped reads

    // you can use this if you don't understand what you are doing...
    bool AG_WILDCARD = false;

    uint32_t max_mismatches = 6; // max mismatches allowed

    /* number of reads to map at one batch */
    uint32_t batch_size = 10000000;
    const uint32_t max_batch_size = 100000000;

    /* ignore seeds which have more than b candidates */
    uint32_t b = 5000;
    uint32_t top_k = 50; // (paired-end) candidates per end keep
    int frag_range = 1000; // max dist btwn concordantly mapped mates
    int num_of_threads = 1; // number of threads to use

    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "map Illumina BS-seq reads","");
    opt_parse.add_opt("index", 'i', "index file created by makedb command "
                      "(the suffix of the index file should be '.dbindex')",
                      true, index_file);
    opt_parse.add_opt("reads", 'r', "comma-sep list of read files "
                      "for singled-end mapping (expect suffix .fastq or .fq)",
                      false, se_read_files_csv);
    opt_parse.add_opt("reads1", '1', "comma-separated list of read files for "
                      "mate 1 (expect suffix .fastq or .fq)",
                      false, pe_read_files_end1_csv);
    opt_parse.add_opt("reads2", '2', "comma-separated list of read files for "
                      "mate 2 (expect suffix .fastq or .fq)",
                      false, pe_read_files_end2_csv);
    opt_parse.add_opt("output", 'o', "output file names (comma sep)",
                      true, output_files_csv);
    opt_parse.add_opt("mismatch", 'm', "max allowed mismatches", false,
                      max_mismatches);
    opt_parse.add_opt("number", 'N', "number of reads per batch",
                      false, batch_size);
    opt_parse.add_opt("ambiguous", 'a', "output one random location for "
                      "ambiguously mapping reads in separate file",
                      false, ambiguous);
    opt_parse.add_opt("unmapped", 'u', "output unmapped reads in separate file",
                      false, unmapped);
    opt_parse.add_opt("clip", 'C', "clip the specified adaptor", false, adaptor);
    opt_parse.add_opt("ag-wild", 'A', "map using A/G bisulfite wildcards "
                      "(single-end)", false, AG_WILDCARD);
    opt_parse.add_opt("bucket", 'b', "maximum candidates for a seed",
                      false, top_k);
    opt_parse.add_opt("topk", 'k', "maximum allowed mappings for a read "
                      "(paired-end)", false, top_k);
    opt_parse.add_opt("fraglen", 'L', "max fragment length (paired-end)", false,
                      frag_range);
    opt_parse.add_opt("sam", '\0', "output sam format", false, sam_format);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("thread", 't', "number of threads for mapping", false,
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
    if (!leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    validate_index_file(index_file);

    vector<string> se_read_files;
    split_filenames(se_read_files_csv, se_read_files);
    for (auto && i : se_read_files)
      if (!valid_suffix(i, fastq_suffixes))
        throw runtime_error("read file invalid suffix: " + i);

    vector<string> pe_read_files_end1;
    split_filenames(pe_read_files_end1_csv, pe_read_files_end1);

    vector<string> pe_read_files_end2;
    split_filenames(pe_read_files_end2_csv, pe_read_files_end2);

    if (pe_read_files_end1.size() != pe_read_files_end2.size())
      throw runtime_error("unequal number of end1 and end2 files");

    for (auto && i : pe_read_files_end1)
      if (!valid_suffix(i, fastq_suffixes))
        throw runtime_error("read file invalid suffix: " + i);

    for (auto && i : pe_read_files_end2)
      if (!valid_suffix(i, fastq_suffixes))
        throw runtime_error("read file invalid suffix: " + i);

    vector<string> output_files;
    split_filenames(output_files_csv, output_files);
    if (output_files.size() !=
        se_read_files.size() + pe_read_files_end1.size())
      throw runtime_error("wrong number of output files: " +
                          output_files_csv);

    //////////////////////////////////////////////////////////////
    // CHECK OPTIONS
    if (VERBOSE)
      cerr << "max_mismatches: " << max_mismatches << endl
           << "threads: " << num_of_threads << endl;

    if (batch_size > max_batch_size)
      throw runtime_error("batch size may not exceed" +
                          to_string(max_batch_size));

    if (top_k < 2 || top_k > 300)
      throw runtime_error("paired-end candidates must be in [2, 300]");

    if (VERBOSE)
      ShowGenomeInfo(index_file);

    // map any single-end reads
    const size_t n_se_read_files = se_read_files.size();
    if (VERBOSE)
      cerr << "n_se_read_files: " << n_se_read_files << endl;
    for (uint32_t i = 0; i < se_read_files.size(); ++i)
      ProcessSingledEndReads(VERBOSE, index_file, se_read_files[i],
                             output_files[i], batch_size,
                             max_mismatches, b, adaptor, AG_WILDCARD,
                             ambiguous, unmapped, sam_format, num_of_threads);

    // now map any paired-end reads
    const size_t n_pe_read_files = pe_read_files_end1.size();
    if (VERBOSE)
      cerr << "n_pe_read_files: " << n_pe_read_files << endl;
    for (uint32_t i = 0; i < n_pe_read_files; ++i)
      ProcessPairedEndReads(VERBOSE, index_file,
                            pe_read_files_end1[i], pe_read_files_end2[i],
                            output_files[i+n_se_read_files],
                            batch_size, max_mismatches, b, adaptor,
                            top_k, frag_range, ambiguous, unmapped,
                            sam_format, num_of_threads);
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
