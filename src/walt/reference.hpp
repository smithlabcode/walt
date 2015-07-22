/*
 *    The is the header file for building index.
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

/*
 * reference.cpp loads chromosomes from the genome file, builds hash table
 * and stores them in struct Genome and HashTable.
 *
 * All the characters in the genome are converted into capital letters.
 *
 * Ns in the genome are converted to T when mapping _1 read files, and
 * converted to A when mapping _2 read files.
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "util.hpp"

#include <set>
#include <vector>
#include <string>
#include <utility>

using std::set;
using std::pair;
using std::string;
using std::vector;
using std::make_pair;

struct Genome {
  /* chromosome name */
  vector<string> name;

  /* chromosome length */
  vector<uint32_t> length;

  /* all chromosomes are concatenated to one string, and stored
   * in vector<char> sequence. start_index indicates the start position
   * of the chromosome in vector<char> sequence. */
  vector<uint32_t> start_index;

  /* There are two different strands, '+' and '-'. The chromosome read from
   * reference file is '+', and the one got from reverse compliment rule
   * is '-'.  */
  char strand;

  /* number of chromosomes in the genome */
  uint32_t num_of_chroms;

  /* The total number of nucleotides in the genome, which is the sum of
   * all chromosomes' length. */
  uint32_t length_of_genome;

  /* all chromosome are concatenated to one string, and stored in sequence. */
  vector<char> sequence;
};

/* HashTable stores genome positions for each k-mer
 * HashTable is changed from
 * unordered_map<uint32_t, vector<GenomePosition> >
 * to current struct in order to reduce the time for reading index. When using
 * unordered_map, the number of reading disk equals to the number of k-mers.
 * Now all the positions are put in one array index, only read disk once.
 * */
struct HashTable {
  /* counter_size records the size of array counter */
  uint32_t counter_size;

  /* index_size records the size of array index */
  uint32_t index_size;

  /* counter is a indicator array. It recodes the start positions of
   * each k-mer in the index array */
  vector<uint32_t> counter;

  /* index array stores genome positions for each k-mer */
  vector<uint32_t> index;
};

/* find the first index in nums which is larger or equal to pos */
uint32_t getChromID(const vector<uint32_t>& nums, const uint32_t& pos);

/* identify chromosome files and estimate the size of each chromosome */
void IdentifyChromosomes(const string& chrom_file, vector<string>& chrom_files);

/* read chromosomes from disk and store in genome */
void ReadGenome(const vector<string>& chrom_files, Genome& genome);

/* get the reverse complimentary strand of genome */
void ReverseComplementGenome(Genome& genome);

/* Cs in the genome are converted to Ts */
void C2T(vector<char>& sequence);

/* Gs in the genome are converted to As */
void G2A(vector<char>& sequence);

/* count how many k-mers for each hash value (bucket) */
void CountBucketSize(const Genome& genome, HashTable& hash_table,
                     set<uint32_t>& extremal_large_bucket);

/* put genome positions to the corresponding bucket */
void HashToBucket(const Genome& genome, HashTable& hash_table,
                  const set<uint32_t>& extremal_large_bucket);

/* Sort each bucket, if the seed length is more than 12, then use binary search
 * for the rest part of the seed */
void SortHashTableBucket(const Genome& genome, HashTable& hash_table);

#ifdef DEBUG
/* Output the Hash Table to a human readable file for testing */
void TestHashTable(const Genome& genome, const HashTable& hash_table);
#endif

/* After building the hash table for all chromosomes, write them to the disk.
 * Next time when mapping the reads, first using ReadIndex function to read
 * the chromosomes and hash tables */
void WriteIndex(const string& index_file, const Genome& genome,
                const HashTable& hash_table);

/* read the chromosomes and hash tables from the disk */
void ReadIndex(const string& index_file, Genome& genome, HashTable& hash_table);

/* write the head information to disk, including 4 index names, chromosome names
 * and lengths, and also the largest size of index array in HashTable */
void WriteIndexHeadInfo(const string& index_file, const Genome& genome,
                        const uint32_t& size_of_index);

/* read the head information from disk */
void ReadIndexHeadInfo(const string& index_file, Genome& genome,
                       uint32_t& size_of_index);

/* show genome length and the number of chromosomes in the genome */
void ShowGenomeInfo(const string& index_file);

/* show SAM head information */
void SAMHead(const string& index_file, const string& command, FILE * fout);

#endif /* REFERENCE_H_ */
