/*
 * This file reads reference genome and stores them in vector<Chromosome>.
 * All the characters in genome are transfered to capital letters, and Ns
 * in the genome are randomly transfered to A, C, G or T.
 * For mapping bisulfite seqeunces, Cs in the genome are transfered to T.
 * */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "util.hpp"

#include <vector>
#include <string>
#include <tr1/unordered_map>

using std::string;
using std::vector;
using std::tr1::unordered_map;

/* HashTable stores genome positions for each k-mer */
typedef std::tr1::unordered_map<uint32_t, vector<uint32_t> > HashTable;

struct Chromosome {
  string name;
  uint32_t length;
  char strand;

  vector<char> sequence;

  HashTable hash_table;
};

typedef vector<Chromosome> Genome;

void IdentifyChromosomes(const string& chrom_file, vector<string>& chrom_files);
void ReadChromsAndBuildIndex(const vector<string>& chrom_files, Genome* genome);
void WriteIndex(const string& index_file, const Genome& genome);

#endif /* REFERENCE_H_ */
