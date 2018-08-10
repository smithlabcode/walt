/*
 *    The head file for mapping single-end reads
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

#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include "reference.hpp"
#include "smithlab_utils.hpp"

#include <string>
#include <vector>

/* BestMatch records the match position with the minimal number
 * of mismatches. times records how many genome positions with
 * this minimal number of mismatches. If times equals 1, means
 * the read is uniquely matched to the genome with the minimal
 * number of mismatches. */
struct BestMatch {
  BestMatch(const uint32_t& _genome_pos = 0, const uint32_t& _times = 0,
            const char& _strand = '+', const uint32_t& _mismatch = MAX_UINT32)
      : genome_pos(_genome_pos),
        times(_times),
        strand(_strand),
        mismatch(_mismatch) {
  }

  uint32_t genome_pos;
  uint32_t times;
  char strand;
  uint32_t mismatch;
};

/* count the number of uniquely mapped, ambiguous mapped and unmapped reads */
struct StatSingleReads {
  StatSingleReads(const bool& _ambiguous, const bool& _unmapped,
                  const std::string &output_file, const bool& _SAM)
      : ambiguous(_ambiguous),
        unmapped(_unmapped),
        SAM(_SAM) {
    total_reads = 0;
    unique_mapped_reads = 0;
    ambiguous_mapped_reads = 0;
    unmapped_reads = 0;

    num_of_short_reads = 0;

    if (ambiguous && !SAM) {
      fambiguous = fopen(std::string(output_file + "_ambiguous").c_str(), "w");
      if (!fambiguous) {
        throw SMITHLABException(
            "cannot open input file " + std::string(output_file + "_ambiguous"));
      }
    }
    if (unmapped && !SAM) {
      funmapped = fopen(std::string(output_file + "_unmapped").c_str(), "w");
      if (!funmapped) {
        throw SMITHLABException(
            "cannot open input file " + std::string(output_file + "_unmapped"));
      }
    }
  }
  ~StatSingleReads() {
    if (ambiguous && !SAM) {
      fclose(fambiguous);
    }
    if (unmapped && !SAM) {
      fclose(funmapped);
    }
  }

  uint32_t total_reads;
  uint32_t unique_mapped_reads;
  uint32_t ambiguous_mapped_reads;
  uint32_t unmapped_reads;

  // number of reads that are shorter than the minimum length requirement
  uint32_t num_of_short_reads;

  FILE * fambiguous;
  FILE * funmapped;

  bool ambiguous;
  bool unmapped;
  bool SAM;
};

/* load reads from reads file, each time load n_reads_per_batch reads,
 * start from  read_start_idx */
void LoadReadsFromFastqFile(FILE * fin, const uint32_t read_start_idx,
                            const uint32_t& n_reads_per_batch,
                            const std::string &adaptor, uint32_t& num_of_reads,
                            std::vector<std::string> &read_names,
                            std::vector<std::string> &read_seqs,
                            std::vector<std::string> &read_scores);

/* reverse the std::string */
std::string ReverseString(const std::string &str);

/* reverse compliment std::string */
std::string ReverseComplimentString(const std::string &str);

/* reads from _1 file, Cs are converted to Ts*/
void C2T(const std::string &org_read, const uint32_t& read_len, std::string &read);

/* reads from _2 file, Gs are converted to As*/
void G2A(const std::string &org_read, const uint32_t& read_len, std::string &read);

/* find the region of index where those positions started with the seed */
void IndexRegion(const std::string &read, const Genome& genome,
                 const HashTable& hash_table, const uint32_t& seed_len,
                 std::pair<uint32_t, uint32_t> &region);

/* output the uniquely mapped reads or ambiguously mapped reads */
void OutputUniquelyAndAmbiguousMapped(const BestMatch& best_match,
                                      const std::string &read_name,
                                      const std::string &read_seq,
                                      const std::string &read_score,
                                      const Genome& genome,
                                      const bool& AG_WILDCARD, FILE * fout);

/* output the unmapped reads */
void OutputUnmapped(const std::string &read_name, const std::string &read_seq,
                    const std::string &read_score, FILE * fout);

/* output the single end results */
void OutputSingleResults(const BestMatch& best_match, const std::string &read_name,
                         const std::string &read_seq, const std::string &read_score,
                         const Genome& genome, const bool& AG_WILDCARD,
                         StatSingleReads& stat_single_reads, FILE * fout);

/* update number of unmapped, uniquely mapped and ambiguously mapped reads */
void StatInfoUpdate(const uint32_t times, StatSingleReads& stat_single_reads);

/* singled-end read */
void ProcessSingledEndReads(const bool VERBOSE, const std::string &index_file,
                            const std::string &input_file,
                            const std::string &output_file,
                            const uint32_t n_reads_per_batch,
                            const uint32_t max_mismatches, const uint32_t b,
                            const std::string &adaptor,
                            const bool AG_WILDCARD, const bool ambiguous,
                            const bool unmapped, const bool SAM,
                            const int num_of_threads);

#endif /* MAPPING_HPP_ */
