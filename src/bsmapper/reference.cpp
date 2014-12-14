/* The detail description of each function please refer to head file */

#include <cmath>

#include "reference.hpp"
#include "smithlab_os.hpp"

#include <fstream>
#include <algorithm>

void IdentifyChromosomes(const string& chrom_file,
                         vector<string>& chrom_files) {
  cerr << "[IDENTIFYING CHROMS] ";
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, "fa", chrom_files);
  } else {
    chrom_files.push_back(chrom_file);
  }

  cerr << "[DONE]" << endl << "chromosome files found (approx size):" << endl;
  for (uint32_t i = 0; i < chrom_files.size(); ++i) {
    cerr << chrom_files[i] << " ("
         << roundf(get_filesize(chrom_files[i]) / 1e06) << "Mbp)" << endl;
  }
  cerr << endl;
}

void CountBucketSize(const Genome& genome, HashTable* hash_table) {
  cerr << "[COUNT HASH BUCKET SIZE]" << endl;
  hash_table->counter_size = 1 << (2 * F2SEEDWIGTH);
  cerr << "[THE SIZE OF COUNTER ARRAY IS " << hash_table->counter_size << "]"
       << endl;
  hash_table->counter.resize(hash_table->counter_size + 1);
  for (uint32_t i = 0; i <= hash_table->counter_size; ++i) {
    hash_table->counter[i] = 0;
  }

  uint32_t size = 0, hash_value = 0;
  for (uint32_t i = 0; i < genome.size(); ++i) {
    if (genome[i].length < HASHLEN)
      continue;
    size = genome[i].length - HASHLEN;
    for (uint32_t j = 0; j <= size; ++j) {
      hash_value = getHashValue(&(genome[i].sequence[j]));
      hash_table->counter[hash_value]++;
    }
  }

  for (uint32_t i = 1; i <= hash_table->counter_size; ++i) {
    hash_table->counter[i] += hash_table->counter[i - 1];
  }
  hash_table->index_size = hash_table->counter[hash_table->counter_size];

  for (uint32_t i = hash_table->counter_size - 1; i >= 1; --i) {
    hash_table->counter[i] = hash_table->counter[i - 1];
  }
  hash_table->counter[0] = 0;
}

void HashToBucket(const Genome& genome, HashTable* hash_table) {
  cerr << "[HASH TO BUCKET]" << endl;
  cerr << "[THE MEMORY OF INDEX ARRAY IS "
       << sizeof(uint32_t) * (hash_table->index_size / GB) << " GB]" << endl;
  hash_table->index.resize(hash_table->index_size);

  uint32_t size = 0, hash_value = 0;
  for (uint32_t i = 0; i < genome.size(); ++i) {
    if (genome[i].length < HASHLEN)
      continue;
    size = genome[i].length - HASHLEN;
    for (uint32_t j = 0; j <= size; ++j) {
      hash_value = getHashValue(&(genome[i].sequence[j]));
      hash_table->index[hash_table->counter[hash_value]++] = GenomePosition(i,
                                                                            j);
    }
  }

  for (uint32_t i = hash_table->counter_size - 1; i >= 1; --i) {
    hash_table->counter[i] = hash_table->counter[i - 1];
  }
  hash_table->counter[0] = 0;
}

struct SortHashTableBucketCMP {
  SortHashTableBucketCMP(const Genome* _genome)
      : genome(_genome) {
  }
  bool operator()(const GenomePosition& p1, const GenomePosition& p2) {
    const char* c_seq1 = &((*genome)[p1.chrom_id].sequence[p1.chrom_pos]);
    const char* c_seq2 = &((*genome)[p2.chrom_id].sequence[p2.chrom_pos]);
    uint32_t l1 = (*genome)[p1.chrom_id].length;
    uint32_t l2 = (*genome)[p2.chrom_id].length;

    for (uint32_t j = F2SEEDWIGTH; j < 32; ++j) {
      if (F2SEEDPAOSITION[j] >= l1)
        return true;
      if (F2SEEDPAOSITION[j] >= l2)
        return false;

      if (c_seq1[F2SEEDPAOSITION[j]] < c_seq2[F2SEEDPAOSITION[j]])
        return true;
      else if (c_seq1[F2SEEDPAOSITION[j]] > c_seq2[F2SEEDPAOSITION[j]])
        return false;
    }
    return false;
  }
  const Genome* genome;
};

/* Sort each bucket, if the seed lenght is more than 12, then use binary search for
 * the left part of the seed */
void SortHashTableBucket(const Genome* genome, HashTable * hash_table) {
  cerr << "[SORTING BUCKETS FOR HASH TABLE] " << endl;
  for (uint32_t i = 0; i < hash_table->counter_size; ++i) {
    std::sort(hash_table->index.begin() + hash_table->counter[i],
              hash_table->index.begin() + hash_table->counter[i + 1],
              SortHashTableBucketCMP(genome));
  }
}

void TestHashTable(const Genome& genome, const HashTable& hash_table) {
  cerr << "[TEST HASH TABLE] " << endl;
  std::ofstream fout("test.txt");
  for (uint32_t i = 0; i <= hash_table.counter_size; ++i) {
    for (uint64_t j = hash_table.counter[i]; j < hash_table.counter[i + 1];
        ++j) {
      GenomePosition id = hash_table.index[j];
      const char* seq = &(genome[id.chrom_id].sequence[id.chrom_pos]);
      for (uint32_t k = 0; k < 32; ++k) {
        fout << seq[F2SEEDPAOSITION[k]];
      }
      fout << " " << i << " " << id.chrom_pos << std::endl;
    }
    fout << "-----------------------------------" << endl;
  }
  fout.close();
}

void BuildHashTable(const Genome& genome, HashTable* hash_table) {
  cerr << "[BUILD HASH TABLE]" << endl;
  TIME_INFO(CountBucketSize(genome, hash_table), "COUNT BUCKET SIZE");
  TIME_INFO(HashToBucket(genome, hash_table), "HASH TO BUCKET");
  TIME_INFO(SortHashTableBucket(&genome, hash_table), "SORT BUCKETS");
}

void ToUpper(Chromosome* chrom) {
  for (uint32_t i = 0; i < chrom->length; ++i) {
    chrom->sequence[i] = toupper(chrom->sequence[i]);
  }
}

void N2ACGT(Chromosome* chrom) {
  srand(time(NULL));
  for (uint32_t i = 0; i < chrom->length; ++i) {
    if ('N' == chrom->sequence[i]) {
      //int r = rand() % 4;
      int r = 3;
      chrom->sequence[i] = getNT(r);
    }
  }
}

void C2T(Chromosome* chrom) {
  for (uint32_t i = 0; i < chrom->length; ++i) {
    if ('C' == chrom->sequence[i]) {
      chrom->sequence[i] = 'T';
    }
  }
}

void ReadGenome(const vector<string>& chrom_files, Genome* genome) {
  cerr << "[READING CHROMOSOMES] " << endl;
  vector<string> chrom_names;
  vector<string> chrom_seqs;
  uint64_t all_chroms_len = 0;
  for (uint32_t i = 0; i < chrom_files.size(); ++i) {
    vector<string> tmp_chrom_names;
    vector<string> tmp_chrom_seqs;

    /* read chromosome name and seqeunce from the chromosome file */
    read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, tmp_chrom_seqs);

    for (uint32_t j = 0; j < tmp_chrom_seqs.size(); ++j) {
      chrom_names.push_back(tmp_chrom_names[j]);
      chrom_seqs.push_back(tmp_chrom_seqs[j]);
      all_chroms_len += tmp_chrom_seqs[j].size();
    }
  }

  /* copy chroms sequences to genome */
  uint32_t num_of_chroms = chrom_seqs.size();
  genome->resize(2 * num_of_chroms);
  cerr << "[THERE ARE " << num_of_chroms << " CHROMOSOMES]" << endl;
  cerr << "[THE TAOTAL LENGTH OF ALL CHROMOSOMES IS " << all_chroms_len << "]"
       << endl;
  cerr << "[USING FIRST " << HASHLEN << " NUCLEOTIDES AS THE HASH KEY]" << endl;
  for (uint32_t i = 0; i < num_of_chroms; ++i) {
    Chromosome& chrom = (*genome)[2 * i];
    Chromosome& chrom_rc = (*genome)[2 * i + 1];

    chrom.name = chrom_names[i];
    chrom.length = chrom_seqs[i].size();
    chrom.strand = '+';
    chrom.sequence.resize(chrom.length);
    for (uint32_t j = 0; j < chrom.length; ++j) {
      chrom.sequence[j] = chrom_seqs[i][j];
    }
    ToUpper(&chrom);
    N2ACGT(&chrom);

    /* reverse compliment strand */
    chrom_rc.name = chrom_names[i];
    chrom_rc.length = chrom_seqs[i].size();
    chrom_rc.strand = '-';
    chrom_rc.sequence.resize(chrom_rc.length);
    for (uint32_t j = 0; j < chrom_rc.length; ++j) {
      chrom_rc.sequence[j] = complimentBase(
          chrom.sequence[chrom.length - j - 1]);
    }
    ToUpper(&chrom_rc);
    N2ACGT(&chrom_rc);

    C2T(&chrom);
    C2T(&chrom_rc);
  }
}

void WriteIndex(const string& index_file, const Genome& genome,
                const HashTable& hash_table) {
  cerr << "[WRITTING INDEX TO " << index_file << "]" << endl;
  FILE * fout = fopen(index_file.c_str(), "wb");

  uint32_t num_of_chroms = genome.size();
  fwrite(&num_of_chroms, sizeof(uint32_t), 1, fout);
  char chrom_name[256];
  for (uint32_t i = 0; i < num_of_chroms; ++i) {
    /* write chromosome to disk */
    uint32_t chrom_name_len = genome[i].name.size();
    if (chrom_name_len > 255) {
      chrom_name_len = 255;
    }
    fwrite(&chrom_name_len, sizeof(uint32_t), 1, fout);
    strcpy(chrom_name, genome[i].name.substr(0, chrom_name_len).c_str());
    fwrite(chrom_name, sizeof(char), chrom_name_len, fout);
    fwrite(&(genome[i].length), sizeof(uint32_t), 1, fout);
    fwrite(&(genome[i].strand), sizeof(char), 1, fout);
    fwrite(&(genome[i].sequence[0]), sizeof(char), genome[i].length, fout);
  }

  /* write hash table to disk */
  fwrite(&(hash_table.counter_size), sizeof(uint32_t), 1, fout);
  fwrite(&(hash_table.index_size), sizeof(uint64_t), 1, fout);
  fwrite(&(hash_table.counter[0]), sizeof(uint64_t),
         hash_table.counter_size + 1, fout);
  fwrite(&(hash_table.index[0]), sizeof(GenomePosition), hash_table.index_size,
         fout);

  fclose(fout);
}

void ReadIndex(const string& index_file, Genome* genome,
               HashTable* hash_table) {
  cerr << "[READING INDEX FROM " << index_file << "]" << endl;
  FILE * fin = fopen(index_file.c_str(), "rb");
  FILE_OPEN_CHECK(fin);

  uint32_t num_of_chroms;
  FREAD_CHECK(fread(&num_of_chroms, sizeof(uint32_t), 1, fin), 1);
  genome->resize(num_of_chroms);
  cerr << "[THERE ARE " << num_of_chroms << " CHROMOSOMES IN THE GENOME]"
       << endl;

  char chrom_strand;
  char chrom_name[256];
  uint32_t chrom_name_len, chrom_length;
  for (uint32_t i = 0; i < num_of_chroms; ++i) {
    cerr << "[" << i + 1 << "/" << num_of_chroms << "]";
    /* read chromosome from disk */
    FREAD_CHECK(fread(&chrom_name_len, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(chrom_name, sizeof(char), chrom_name_len, fin),
                chrom_name_len);
    chrom_name[chrom_name_len] = 0;
    FREAD_CHECK(fread(&chrom_length, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(&chrom_strand, sizeof(char), 1, fin), 1);
    (*genome)[i].name = chrom_name;
    (*genome)[i].length = chrom_length;
    (*genome)[i].strand = chrom_strand;

    (*genome)[i].sequence.resize(chrom_length);
    FREAD_CHECK(
        fread(&((*genome)[i].sequence[0]), sizeof(char), chrom_length, fin),
        chrom_length);
  }
  cerr << endl;

  /* read hash table from disk */

  cerr << "[READING HASH TABLE] ";
  FREAD_CHECK(fread(&(hash_table->counter_size), sizeof(uint32_t), 1, fin), 1);
  FREAD_CHECK(fread(&(hash_table->index_size), sizeof(uint64_t), 1, fin), 1);

  hash_table->counter.resize(hash_table->counter_size + 1);
  hash_table->index.resize(hash_table->index_size);

  FREAD_CHECK(
      fread(&(hash_table->counter[0]), sizeof(uint64_t),
            hash_table->counter_size + 1, fin),
      hash_table->counter_size + 1);
  FREAD_CHECK(
      fread(&(hash_table->index[0]), sizeof(GenomePosition),
            hash_table->index_size, fin),
      hash_table->index_size);

  fclose(fin);
}
