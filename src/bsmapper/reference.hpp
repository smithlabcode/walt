/*
 * This file reads reference genome, builds hash table and stores them
 * in vector<Chromosome> which is typedef as Genome.
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

struct Chromosome {
  /* chromoseome name */
  string name;

  /* chromoseome length */
  uint32_t length;

  /* There are two different strands, '+' and '-'. The one read from reference
   * file is '+', and the one got from reverse compliment rule is '-'.  */
  char strand;

  /* chromosome sequence */
  vector<char> sequence;
};

/* Genome contains several Chromosomes */
typedef vector<Chromosome> Genome;

struct GenomePosition {
  GenomePosition(const uint32_t& _chrom_id, const uint32_t& _chrom_pos)
      : chrom_id(_chrom_id),
        chrom_pos(_chrom_pos) {
  }
  GenomePosition() {
    chrom_id = 0;
    chrom_pos = 0;
  }
  uint32_t chrom_id;
  uint32_t chrom_pos;
};

/*
 * HashTable stores positions in the genome for each k-mer
 */
struct HashTable {
  /* counter_size is the size of counter array */
  uint32_t counter_size;

  /* index_size is the size of index array */
  uint64_t index_size;

  /* counter is an array which indicates the start position of k-mers in index array.
   * So the hash key is the k-mer (transfered to integer), and the hash value is the
   * start position for the partiular k-mer in the index array */
  vector<uint64_t> counter;

  /* index array stores positions for all k-mer */
  vector<GenomePosition> index;
};

/* identify all the chromosome files and estimate the size of each chromosome */
void IdentifyChromosomes(const string& chrom_file, vector<string>& chrom_files);

/* get the data for the struct of Chromosome, first read the chromosome sequence
 * from the chrom_files, then build hash tables for them, each chromosome has their
 * own hash table */
void ReadGenome(const vector<string>& chrom_files, Genome* genome);

/* build hash table contains three parts, first count how many position for each
 * k-mer, and then put the position to the index array, and finally sort each
 * bucket */
void BuildHashTable(const Genome& genome, HashTable* hash_table);

/* Output the Hash Table to a human readable file for testing */
void TestHashTable(const Genome& genome, const HashTable& hash_table);

/* After building the hash table for all the chromosomes, write them to the disk.
 * Next time when mapping the reads, first should using ReadIndex function to read
 * the chromosomes and hash tables */
void WriteIndex(const string& index_file, const Genome& genome,
                const HashTable& hash_table);

/* read the chromosomes and hash tables from the disk */
void ReadIndex(const string& index_file, Genome* genome, HashTable* hash_table);

#endif /* REFERENCE_H_ */
