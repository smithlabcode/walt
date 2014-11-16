/*
 * This file builds index (hash table) for the reference genome.
 * */

#ifndef INDEX_H_
#define INDEX_H_

#include "util.hpp"
#include "reference.hpp"

#include <utility>
#include <string>
#include <tr1/unordered_map>

using std::pair;
using std::string;
using std::make_pair;
using std::tr1::unordered_multimap;

/* HashTable stores positions in the genome for each k-mer */
typedef std::tr1::unordered_map<uint32_t, vector<uint32_t> > HashTable;

class BuildIndex {
 public:
  BuildIndex(const string& _index_file, const Genome* _genome,
             HashTable* _hash_table)
      : index_file(_index_file),
      genome(_genome),
      hash_table(_hash_table) {
        TIME_INFO(BuildHashTable(), "BUILD HASH TABLE");
        WriteIndex();
      }

 private:
  /* build the hash table, keys are k-mers, and values
   * are the positions in the genome */
  void BuildHashTable();

  /* write index to disk */
  void WriteIndex();

  string index_file;
  const Genome* genome;
  HashTable* hash_table;
};

class ReadIndex {
 public:
  ReadIndex(const string& _index_file, Genome* _genome, HashTable* _hash_table);

 private:
  string index_file;
  Genome* genome;
  HashTable* hash_table;
};

#endif /* INDEX_H_ */
