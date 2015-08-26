/*
 *    This is the main function to mark repeats in each hash bucket.
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

struct SortByNucleotides {
  SortByNucleotides(const vector<Genome>& _genome,
                    const uint32_t& _num_of_nucleotides)
      : genome(_genome),
        num_of_nucleotides(_num_of_nucleotides) {
  }
  bool operator()(const pair<uint32_t, int>& p1,
                  const pair<uint32_t, int>& p2) {
    const char* c_seq1 = &(genome[p1.second].sequence[p1.first]);
    const char* c_seq2 = &(genome[p2.second].sequence[p2.first]);

    uint32_t chr_id1 = getChromID(genome[p1.second].start_index, p1.first);
    uint32_t chr_id2 = getChromID(genome[p2.second].start_index, p2.first);

    uint32_t l1 = genome[p1.second].start_index[chr_id1 + 1] - p1.first;
    uint32_t l2 = genome[p2.second].start_index[chr_id2 + 1] - p2.first;

    for (uint32_t j = 0; j < num_of_nucleotides; ++j) {
      /*Strict Weak Ordering */
      if (j >= l2)
        return false;
      if (j >= l1)
        return true;

      if (c_seq1[j] < c_seq2[j])
        return true;
      else if (c_seq1[j] > c_seq2[j])
        return false;
    }
    return false;
  }

  const vector<Genome>& genome;
  const uint32_t& num_of_nucleotides;
};

bool IsRepeats(const vector<Genome>& genome,
               const vector<pair<uint32_t, int> >& bucket,
               const uint32_t& bucket_size, const uint32_t& num_of_nucleotides,
               const uint32_t& pos) {
  if (pos + 1 >= bucket_size)
    return false;
  pair<uint32_t, int> p1 = bucket[pos];
  pair<uint32_t, int> p2 = bucket[pos + 1];

  const char* c_seq1 = &(genome[p1.second].sequence[p1.first]);
  const char* c_seq2 = &(genome[p2.second].sequence[p2.first]);

  uint32_t chr_id1 = getChromID(genome[p1.second].start_index, p1.first);
  uint32_t chr_id2 = getChromID(genome[p2.second].start_index, p2.first);

  uint32_t l1 = genome[p1.second].start_index[chr_id1 + 1] - p1.first;
  uint32_t l2 = genome[p2.second].start_index[chr_id2 + 1] - p2.first;
  if (l1 < num_of_nucleotides || l2 < num_of_nucleotides)
    return false;
  for (uint32_t i = 0; i < num_of_nucleotides; ++i) {
    if (c_seq1[i] != c_seq2[i])
      return false;
  }

  return true;
}

// set bit to 1
void SetBit(const pair<uint32_t, int>& pos, vector<vector<uint32_t> >& mark) {
  uint32_t word = pos.first / 32;
  uint32_t bit = pos.first % 32;

  mark[pos.second][word] |= 1 << bit;
}

void MarkRepeats(const string& index_file, const vector<Genome>& genome,
                 vector<HashTable>& hash_table,
                 const uint32_t& num_of_nucleotide_for_sort) {
  uint32_t num_of_bucket = power(4, F2SEEDWIGTH);
  uint32_t max_bucket_size = 0;
  max_bucket_size =
      hash_table[0].counter[0] + hash_table[1].counter[0] > max_bucket_size ?
          hash_table[0].counter[0] + hash_table[1].counter[0] : max_bucket_size;

  for (uint32_t i = 1; i <= num_of_bucket; ++i) {
    uint32_t bucket_size = 0;
    bucket_size += hash_table[0].counter[i] - hash_table[0].counter[i - 1];
    bucket_size += hash_table[1].counter[i] - hash_table[1].counter[i - 1];
    max_bucket_size =
        bucket_size > max_bucket_size ? bucket_size : max_bucket_size;
  }

  // MARK REPEATS
  uint32_t mark_size = genome[0].length_of_genome / 32 + 1;
  vector<vector<uint32_t> > mark_90(2, vector<uint32_t>(mark_size, 0));
  vector<vector<uint32_t> > mark_100(2, vector<uint32_t>(mark_size, 0));

  uint32_t bucket_size = 0;
  vector<pair<uint32_t, int> > bucket(max_bucket_size);
  for (uint32_t i = 0; i < num_of_bucket; ++i) {
    bucket_size = 0;
    for (int pi = 0; pi < 2; ++pi) {
      for (uint32_t j = hash_table[pi].counter[i];
          j < hash_table[pi].counter[i + 1]; ++j) {
        bucket[bucket_size++] = make_pair(hash_table[pi].index[j], pi);
      }
    }
    std::sort(bucket.begin(), bucket.begin() + bucket_size,
              SortByNucleotides(genome, num_of_nucleotide_for_sort));
    for (uint32_t j = 0; j < bucket_size; ++j) {
      if (IsRepeats(genome, bucket, bucket_size, 90, j)) {
        SetBit(bucket[j], mark_90);
        SetBit(bucket[j + 1], mark_90);
      }
      if (IsRepeats(genome, bucket, bucket_size, 100, j)) {
        SetBit(bucket[j], mark_100);
        SetBit(bucket[j + 1], mark_100);
      }
    }
  }

  FILE * fout = fopen(string(index_file + ".markrp_90").c_str(), "wb");
  if (!fout) {
    throw SMITHLABException(
        "cannot open input file " + string(index_file + ".markrp_90"));
  }
  fwrite(&bucket_size, sizeof(uint32_t), 1, fout);
  fwrite(&(mark_90[0]), sizeof(uint32_t), bucket_size, fout);
  fwrite(&(mark_90[1]), sizeof(uint32_t), bucket_size, fout);
  fclose(fout);

  fout = fopen(string(index_file + ".markrp_100").c_str(), "wb");
  if (!fout) {
    throw SMITHLABException(
        "cannot open input file " + string(index_file + ".markrp_100"));
  }
  fwrite(&bucket_size, sizeof(uint32_t), 1, fout);
  fwrite(&(mark_100[0]), sizeof(uint32_t), bucket_size, fout);
  fwrite(&(mark_100[1]), sizeof(uint32_t), bucket_size, fout);
  fclose(fout);
}

int main(int argc, const char **argv) {
  try {
    /* index file*/
    string index_file;

    /* maximal number of nucleotide to compare as repeats */
    uint32_t num_of_nucleotide_for_sort = 100;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "mark repeated subsequences in each hash bucket",
                           "");
    opt_parse.add_opt(
        "index",
        'i',
        "index file created by makedb command \
        (the suffix of the index file should be '.dbindex')",
        true, index_file);
    opt_parse.add_opt("len", 'l',
                      "maximal number of nucleotides to compare as repeats')",
                      false, num_of_nucleotide_for_sort);

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
    /****************** END COMMAND LINE OPTIONS *****************/

    // LOAD THE INDEX HEAD INFO
    vector<Genome> genome(2);
    vector<HashTable> hash_table(2);

    uint32_t size_of_index;
    for (int i = 0; i < 2; ++i) {
      ReadIndexHeadInfo(index_file, genome[i], size_of_index);
      genome[i].sequence.resize(genome[i].length_of_genome);
      hash_table[i].counter.resize(power(4, F2SEEDWIGTH) + 1);
      hash_table[i].index.resize(size_of_index);
    }

    ////////// MARK REPEATS FOR (C->T)
    ReadIndex(index_file + "_CT00", genome[0], hash_table[0]);
    ReadIndex(index_file + "_CT01", genome[1], hash_table[1]);
    MarkRepeats(index_file, genome, hash_table, num_of_nucleotide_for_sort);

    ////////// MARK REPEATS FOR (G->A)
    ReadIndex(index_file + "_GA10", genome[0], hash_table[0]);
    ReadIndex(index_file + "_GA11", genome[1], hash_table[1]);
    MarkRepeats(index_file, genome, hash_table, num_of_nucleotide_for_sort);
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
