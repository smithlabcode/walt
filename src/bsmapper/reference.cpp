/* The detail description of each function please refer to head file */

#include <cmath>

#include "reference.hpp"
#include "smithlab_os.hpp"

#include <fstream>
#include <algorithm>

/* find the index in array numbs which is first larger or equal to pos */
uint32_t getChromID(const vector<uint32_t>& nums, const uint32_t& pos) {
  uint32_t size = nums.size();
  if (size == 0) {
    fprintf(stderr, "Please Check the start positions of Chromosomes.\n");
    exit(EXIT_FAILURE);
  }

  uint32_t l = 0, h = size - 1;
  while (l < h) {
    uint32_t m = (l + h + 1) / 2;
    if (pos >= nums[m])
      l = m;
    else
      h = m - 1;
  }

  return l;
}

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

void ReadGenome(const vector<string>& chrom_files, Genome* genome) {
  cerr << "[READING CHROMOSOMES] " << endl;
  vector<string> chrom_names;
  vector<string> chrom_seqs;
  uint32_t all_chroms_len = 0;
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
  genome->num_of_chroms = chrom_seqs.size();
  genome->length_of_genome = all_chroms_len;
  cerr << "[THERE ARE " << genome->num_of_chroms << " CHROMOSOMES IN THE GENOME]" << endl;
  cerr << "[THE TAOTAL LENGTH OF ALL CHROMOSOMES IS " << all_chroms_len << "]"
       << endl;
  genome->name.resize(genome->num_of_chroms);
  genome->length.resize(genome->num_of_chroms);

  /* The last element of genome->start-index is the length of all chroms */
  genome->start_index.resize(genome->num_of_chroms + 1);
  genome->strand = '+';
  genome->sequence.resize(all_chroms_len);
  uint32_t k = 0;
  for (uint32_t i = 0; i < genome->num_of_chroms; ++i) {
    genome->name[i] = chrom_names[i];
    genome->length[i] = chrom_seqs[i].size();
    genome->start_index[i] = k;

    for (uint32_t j = 0; j < chrom_seqs[i].size(); ++j) {
      genome->sequence[k++] = toupper(chrom_seqs[i][j]);
    }
  }
  genome->start_index[genome->num_of_chroms] = k;
}

void C2T(vector<char>& sequence) {
  for (uint32_t i = 0; i < sequence.size(); ++i) {
    if ('C' == sequence[i] || 'N' == sequence[i]) {
      sequence[i] = 'T';
    }
  }
}

void A2G(vector<char>& sequence) {
  for (uint32_t i = 0; i < sequence.size(); ++i) {
    if ('A' == sequence[i] || 'N' == sequence[i]) {
      sequence[i] = 'G';
    }
  }
}

void ReverseGenome(const Genome& genome, Genome* rc_genome) {
  rc_genome->num_of_chroms = genome.num_of_chroms;
  rc_genome->name.resize(rc_genome->num_of_chroms);
  rc_genome->length.resize(rc_genome->num_of_chroms);
  rc_genome->start_index.resize(rc_genome->num_of_chroms + 1);
  rc_genome->length_of_genome = genome.length_of_genome;
  rc_genome->strand = '-';
  rc_genome->sequence.resize(rc_genome->length_of_genome);
  for (uint32_t i = 0; i < rc_genome->num_of_chroms; ++i) {
    rc_genome->name[i] = genome.name[i];
    rc_genome->length[i] = genome.length[i];
    rc_genome->start_index[i] = genome.start_index[i];
    for (uint32_t j = 0; j < genome.length[i]; ++j) {
      rc_genome->sequence[j + genome.start_index[i]] = complimentBase(
          genome.sequence[genome.start_index[i + 1] - j - 1]);
    }
  }
  rc_genome->start_index[rc_genome->num_of_chroms] = genome.start_index[genome
      .num_of_chroms];
}

void CountBucketSize(const Genome& genome, HashTable* hash_table) {
  cerr << "[COUNT BUCKET SIZE]" << endl;
  hash_table->counter_size = power(4, F2SEEDWIGTH);
  hash_table->counter.resize(hash_table->counter_size + 1, 0);

  uint32_t size = 0, hash_value = 0;
  for (uint32_t i = 0; i < genome.num_of_chroms; ++i) {
    if (genome.length[i] < HASHLEN)
      continue;
    size = genome.start_index[i + 1] - HASHLEN;
    for (uint32_t j = genome.start_index[i]; j < size; ++j) {
      hash_value = getHashValue(&(genome.sequence[j]));
      hash_table->counter[hash_value]++;
    }
  }

  //////////////////////////////////////////////////////
  // Erase Extremal Large Bucket
  for (uint32_t i = 0; i < hash_table->counter_size; ++i) {
    if (hash_table->counter[i] >= 500000) {
      cerr << "ERASE THE BUCKET " << i << " SINCE ITS SIZE IS "
           << hash_table->counter[i] << endl;
      hash_table->counter[i] = 0;
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
  hash_table->index.resize(hash_table->index_size, 0);

  uint32_t size = 0, hash_value = 0;
  for (uint32_t i = 0; i < genome.num_of_chroms; ++i) {
    if (genome.length[i] < HASHLEN)
      continue;
    size = genome.start_index[i + 1] - HASHLEN;
    for (uint32_t j = genome.start_index[i]; j < size; ++j) {
      hash_value = getHashValue(&(genome.sequence[j]));
      /* Extremal Large Bucket IS DELETED */
      if (hash_table->counter[hash_value]
          == hash_table->counter[hash_value + 1])
        continue;
      hash_table->index[hash_table->counter[hash_value]++] = j;
    }
  }

  for (uint32_t i = hash_table->counter_size - 1; i >= 1; --i) {
    hash_table->counter[i] = hash_table->counter[i - 1];
  }
  hash_table->counter[0] = 0;
}

struct SortHashTableBucketCMP {
  explicit SortHashTableBucketCMP(const Genome* _genome)
      : genome(_genome) {
  }
  bool operator()(const uint32_t& p1, const uint32_t& p2) {
    const char* c_seq1 = &(genome->sequence[p1]);
    const char* c_seq2 = &(genome->sequence[p2]);

    uint32_t chrom_id1 = getChromID(genome->start_index, p1);
    uint32_t chrom_id2 = getChromID(genome->start_index, p2);

    uint32_t l1 = genome->start_index[chrom_id1 + 1] - p1;
    uint32_t l2 = genome->start_index[chrom_id2 + 1] - p2;

    for (uint32_t j = F2SEEDWIGTH; j < F2SEEDPOSITION_SIZE; ++j) {
      if (F2SEEDPOSITION[j] >= l1)
        return true;
      if (F2SEEDPOSITION[j] >= l2)
        return false;

      if (c_seq1[F2SEEDPOSITION[j]] < c_seq2[F2SEEDPOSITION[j]])
        return true;
      else if (c_seq1[F2SEEDPOSITION[j]] > c_seq2[F2SEEDPOSITION[j]])
        return false;
    }
    return false;
  }
  const Genome* genome;
};

void SortHashTableBucket(const Genome* genome, HashTable* hash_table) {
  cerr << "[SORTING BUCKETS FOR HASH TABLE] " << endl;

  for (uint32_t i = 0; i < hash_table->counter_size; ++i) {
    if (hash_table->counter[i + 1] - hash_table->counter[i] <= 1)
      continue;

    std::sort(hash_table->index.begin() + hash_table->counter[i],
              hash_table->index.begin() + hash_table->counter[i + 1],
              SortHashTableBucketCMP(genome));
  }
}

void TestHashTable(const Genome& genome, const HashTable& hash_table) {
  cerr << "[TEST HASH TABLE] " << endl;
  std::ofstream fout("test.txt");
  for (uint32_t i = 0; i < hash_table.counter_size; ++i) {
    for (uint32_t j = hash_table.counter[i]; j < hash_table.counter[i + 1];
        ++j) {
      fout << i << " ";
      const char* seq = &(genome.sequence[hash_table.index[j]]);
      for (uint32_t k = 0; k < 32; ++k) {
        fout << seq[F2SEEDPOSITION[k]];
      }
      fout << " ";
      for (uint32_t k = 0; k < 32; ++k) {
        fout << seq[k];
      }
      fout << " " << j << " " << hash_table.index[j] << std::endl;
    }
    //fout << "-----------------------------------" << endl;
  }
  fout.close();
}

void WriteIndex(const string& index_file, const Genome& genome,
                const HashTable& hash_table) {
  cerr << "[WRITTING INDEX TO " << index_file << "]" << endl;
  FILE * fout = fopen(index_file.c_str(), "wb");

  fwrite(&(genome.strand), sizeof(char), 1, fout);
  fwrite(&(genome.sequence[0]), sizeof(char), genome.length_of_genome, fout);

  /* write hash table to disk */
  fwrite(&(hash_table.counter_size), sizeof(uint32_t), 1, fout);
  fwrite(&(hash_table.index_size), sizeof(uint32_t), 1, fout);

  fwrite(&(hash_table.counter[0]), sizeof(uint32_t),
         hash_table.counter_size + 1, fout);
  fwrite(&(hash_table.index[0]), sizeof(uint32_t), hash_table.index_size, fout);

  fclose(fout);
}

void ReadIndex(const string& index_file, Genome* genome,
               HashTable* hash_table) {
  cerr << "[READING INDEX FROM " << index_file << "]" << endl;
  FILE * fin = fopen(index_file.c_str(), "rb");
  FILE_OPEN_CHECK(fin);

  FREAD_CHECK(fread(&(genome->strand), sizeof(char), 1, fin), 1);
  FREAD_CHECK(
      fread(&(genome->sequence[0]), sizeof(char), genome->length_of_genome,
            fin),
      genome->length_of_genome);

  /* read hash table from disk */
  FREAD_CHECK(fread(&(hash_table->counter_size), sizeof(uint32_t), 1, fin), 1);
  FREAD_CHECK(fread(&(hash_table->index_size), sizeof(uint32_t), 1, fin), 1);

  FREAD_CHECK(
      fread(&(hash_table->counter[0]), sizeof(uint32_t),
            hash_table->counter_size + 1, fin),
      hash_table->counter_size + 1);
  FREAD_CHECK(
      fread(&(hash_table->index[0]), sizeof(uint32_t), hash_table->index_size,
            fin),
      hash_table->index_size);

  fclose(fin);
}

void WriteIndexHeadInfo(const string& index_file,
                        const vector<string>& index_names, const Genome& genome,
                        const uint32_t& size_of_index) {
  FILE * fout = fopen(index_file.c_str(), "wb");

  uint32_t num_of_chroms = genome.num_of_chroms;
  fwrite(&num_of_chroms, sizeof(uint32_t), 1, fout);

  for (uint32_t i = 0; i < num_of_chroms; ++i) {
    uint32_t chrom_name_len = genome.name[i].size();
    if (chrom_name_len > 255) {
      chrom_name_len = 255;
    }
    fwrite(&chrom_name_len, sizeof(uint32_t), 1, fout);
    fwrite(genome.name[i].c_str(), sizeof(char), chrom_name_len, fout);
  }

  fwrite(&(genome.length[0]), sizeof(uint32_t), num_of_chroms, fout);
  fwrite(&(genome.length_of_genome), sizeof(uint32_t), 1, fout);

  uint32_t num_of_index = index_names.size();
  fwrite(&num_of_index, sizeof(uint32_t), 1, fout);
  for (uint32_t i = 0; i < index_names.size(); ++i) {
    uint32_t name_len = index_names[i].size();
    fwrite(&name_len, sizeof(uint32_t), 1, fout);
    fwrite(index_names[i].c_str(), sizeof(char), name_len, fout);
  }

  fwrite(&size_of_index, sizeof(uint32_t), 1, fout);

  fclose(fout);
}

void ReadIndexHeadInfo(const string& index_file, vector<string>* index_names,
                       Genome* genome, uint32_t* size_of_index) {
  FILE * fin = fopen(index_file.c_str(), "rb");
  FILE_OPEN_CHECK(fin);

  uint32_t num_of_chroms;
  FREAD_CHECK(fread(&num_of_chroms, sizeof(uint32_t), 1, fin), 1);
  cerr << "[THERE ARE " << num_of_chroms << " CHROMOSOMES IN THE GENOME]"
       << endl;
  genome->num_of_chroms = num_of_chroms;

  char chrom_name[256];
  uint32_t chrom_name_len;
  for (uint32_t i = 0; i < num_of_chroms; ++i) {
    FREAD_CHECK(fread(&chrom_name_len, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(chrom_name, sizeof(char), chrom_name_len, fin),
                chrom_name_len);
    chrom_name[chrom_name_len] = 0;
    genome->name[i] = chrom_name;
  }
  FREAD_CHECK(fread(&(genome->length[0]), sizeof(uint32_t), num_of_chroms, fin),
              num_of_chroms);
  genome->start_index[0] = 0;
  for (uint32_t i = 1; i <= num_of_chroms; ++i) {
    genome->start_index[i] = genome->start_index[i - 1] + genome->length[i - 1];
  }

  FREAD_CHECK(fread(&(genome->length_of_genome), sizeof(uint32_t), 1, fin), 1);

  uint32_t num_of_index, name_len;
  char index_name[1024];
  FREAD_CHECK(fread(&num_of_index, sizeof(uint32_t), 1, fin), 1);
  index_names->clear();
  for (uint32_t i = 0; i < num_of_index; ++i) {
    FREAD_CHECK(fread(&name_len, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(index_name, sizeof(char), name_len, fin), name_len);
    index_name[name_len] = 0;
    index_names->push_back(chrom_name);
  }

  FREAD_CHECK(fread(size_of_index, sizeof(uint32_t), 1, fin), 1);

  fclose(fin);
}
