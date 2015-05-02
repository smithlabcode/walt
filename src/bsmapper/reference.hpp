/*
 * This file reads reference genome, builds hash table and stores them
 * in Genome and HashTable.
 * All the characters in genome are transfered into capital letters, and
 * Ns in the genome are randomly transfered to A, C, G or T.
 * For mapping bisulfite seqeunces, Cs in the genome (both strands) are
 * transfered to T.
 * */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "util.hpp"

#include <vector>
#include <string>
#include <utility>

using std::string;
using std::vector;
using std::pair;
using std::make_pair;

struct Genome {
  /* chromosome name */
  vector<string> name;

  /* chromosome length */
  vector<uint32_t> length;

  /* all chromosomes are concatenated to one string, and stored
   * in vector<char> sequence. start_index indicates the start position
   * of the chromosome in vector<char> sequence */
  vector<uint32_t> start_index;

  /* There are two different strands, '+' and '-'. The read from reference
   * file is '+', and the one got from reverse compliment rule is '-'.  */
  char strand;

  /* number of chromosomes in the genome */
  uint32_t num_of_chroms;

  /* The total number of nucleotides in the genome, which is the sum of
   * all chromosomes length*/
  uint32_t length_of_genome;

  /* all chromosome are concatenated to one string, and stored in sequence */
  vector<char> sequence;
};

/* HashTable stores genome positions for each k-mer
 * HashTable is changed from
 * unordered_map<uint32_t, vector<GenomePosition> >
 * to current struct in order to reduce the time for reading index. When using
 * unordered_map, the number of reading disk equals to the number of k-mers.
 * Now we put all the positions in one array index, so we only read disk once.
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

/* find the position of the number which is first larger or equal to pso */
uint32_t getChromID(const vector<uint32_t>& nums, const uint32_t& pos);

/* identify all the chromosome files and estimate the size of each chromosome */
void IdentifyChromosomes(const string& chrom_file, vector<string>& chrom_files);

/* get the reverse complimentary strand of genome*/
void ReverseGenome(const Genome& genome, Genome* rc_genome);

/* read chroms from disk and store in genome */
void ReadGenome(const vector<string>& chrom_files, Genome* genome);

/* Cs in read and genome are transferred to Ts */
void C2T(vector<char>& sequence);

/* As in read and genome are transferred to Gs */
void A2G(vector<char>& sequence);

/* Sort each bucket, if the seed lenght is more than 12, then use binary search for
 * the left part of the seed */
void SortHashTableBucket(const Genome* genome, HashTable * hash_table);

/* Output the Hash Table to a human readable file for testing */
void TestHashTable(const Genome& genome, const HashTable& hash_table);

void CountBucketSize(const Genome& genome, HashTable* hash_table);
void HashToBucket(const Genome& genome, HashTable* hash_table);

/* After building the hash table for all the chromosomes, write them to the disk.
 * Next time when mapping the reads, first should using ReadIndex function to read
 * the chromosomes and hash tables */
void WriteIndex(const string& index_file, const Genome& genome,
                const HashTable& hash_table);

/* read the chromosomes and hash tables from the disk */
void ReadIndex(const string& index_file, Genome* genome, HashTable* hash_table);

/* write the head information to disk, including 4 indexes name, genome name
 * and length, and also the largest size of index array in HashTable */
void WriteIndexHeadInfo(const string& index_file, const Genome& genome,
                        const uint32_t& size_of_index);

/* read the head information from disk */
void ReadIndexHeadInfo(const string& index_file, Genome* genome,
                       uint32_t* size_of_index);

#endif /* REFERENCE_H_ */
