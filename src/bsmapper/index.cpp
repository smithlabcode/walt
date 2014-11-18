/* The detail description of each function please refer to head file */
#include <set>

#include "index.hpp"

using std::set;

void BuildIndex::BuildHashTable() {
  cerr << "[BUILD HASH TABLE]" << endl;
  uint32_t size = 0, hashValue = 0;
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    cerr << "[" << i + 1 << "/" << genome->num_of_chroms << "]";
    if (genome->chrom_sizes[i] + genome->chrom_start_pos[i] < HASHLEN)
      continue;
    size = genome->chrom_sizes[i] + genome->chrom_start_pos[i] - HASHLEN;
    for (uint32_t j = genome->chrom_start_pos[i]; j <= size; ++j) {
      hashValue = getHashValue(&(genome->chrom_seqs[j]));
      (*hash_table)[hashValue].push_back(j);
    }
  }
}

void BuildIndex::WriteIndex() {
  FILE * fout = fopen(index_file.c_str(), "wb");
  cerr << "[WRITTING INDEX TO " << index_file << "]" << endl;

  /* write genome to disk */
  uint32_t chrom_name_len;
  fwrite(&genome->num_of_chroms, sizeof(uint16_t), 1, fout);
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    chrom_name_len = genome->chrom_names[i].size();
    if (chrom_name_len > 255)
      chrom_name_len = 255;
    fwrite(&chrom_name_len, sizeof(uint32_t), 1, fout);
    fwrite(&(genome->chrom_names[i][0]), sizeof(char), chrom_name_len, fout);
    fwrite(&(genome->chrom_sizes[i]), sizeof(uint32_t), 1, fout);
    fwrite(&(genome->chrom_start_pos[i]), sizeof(uint32_t), 1, fout);
  }
  fwrite(&(genome->chrom_seqs[0]), sizeof(char), genome->all_chroms_len, fout);

  /* write hash table to disk */
  uint32_t num_of_keys = hash_table->size();
  fwrite(&(num_of_keys), sizeof(uint32_t), 1, fout);
  for (HashTable::const_iterator it = hash_table->begin();
      it != hash_table->end(); ++it) {
    uint32_t hash_key = it->first;
    fwrite(&(hash_key), sizeof(uint32_t), 1, fout);
    uint32_t num_of_values = it->second.size();
    fwrite(&(num_of_values), sizeof(uint32_t), 1, fout);
    fwrite(&(it->second[0]), sizeof(uint32_t), num_of_values, fout);
  }

  fclose(fout);
}

ReadIndex::ReadIndex(const string& _index_file, Genome* _genome,
                     HashTable* _hash_table)
    : index_file(_index_file),
      genome(_genome),
      hash_table(_hash_table) {

  cerr << "[READING INDEX FROM " << index_file << "]" << endl;
  FILE * fin = fopen(index_file.c_str(), "rb");
  FILE_OPEN_CHECK(fin);

  /* read genome from disk */
  FREAD_CHECK(fread(&genome->num_of_chroms, sizeof(uint16_t), 1, fin), 1);
  uint32_t chrom_name_len, chrom_size, chrom_start_pos;
  char chrom_name[256];
  genome->all_chroms_len = 0;
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    FREAD_CHECK(fread(&chrom_name_len, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(chrom_name, sizeof(char), chrom_name_len, fin),
                chrom_name_len);
    chrom_name[chrom_name_len] = 0;
    genome->chrom_names.push_back(chrom_name);
    FREAD_CHECK(fread(&chrom_size, sizeof(uint32_t), 1, fin), 1);
    genome->chrom_sizes.push_back(chrom_size);
    FREAD_CHECK(fread(&chrom_start_pos, sizeof(uint32_t), 1, fin), 1);
    genome->chrom_start_pos.push_back(chrom_start_pos);
    genome->all_chroms_len += genome->chrom_sizes[i];
  }
  genome->chrom_seqs.resize(genome->all_chroms_len);
  FREAD_CHECK(
      fread(&(genome->chrom_seqs[0]), sizeof(char), genome->all_chroms_len,
            fin),
      genome->all_chroms_len);

  /* read hash table from disk */
  uint32_t num_of_keys = 0, hash_key = 0, num_of_values = 0;
  FREAD_CHECK(fread(&num_of_keys, sizeof(uint32_t), 1, fin), 1);
  for (uint32_t i = 0; i < num_of_keys; ++i) {
    FREAD_CHECK(fread(&hash_key, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(&num_of_values, sizeof(uint32_t), 1, fin), 1);
    vector < uint32_t > hash_values(num_of_values);
    FREAD_CHECK(fread(&(hash_values[0]), sizeof(uint32_t), num_of_values, fin),
                num_of_values);

    hash_table->insert(make_pair(hash_key, hash_values));
  }

  fclose(fin);
}
