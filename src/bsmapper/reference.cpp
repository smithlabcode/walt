/* The detail description of each function please refer to head file */

#include <cmath>

#include "reference.hpp"
#include "smithlab_os.hpp"

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

void BuildHashTable(Chromosome* chrom) {
  uint32_t size = 0, hash_value = 0;
  if (chrom->length < HASHLEN)
    return;
  size = chrom->length - HASHLEN;
  for (uint32_t j = 0; j <= size; ++j) {
    hash_value = getHashValue(&(chrom->sequence[j]));
    chrom->hash_table[hash_value].push_back(j);
  }
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
      int r = rand() % 4;
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

void ReadChromsAndBuildIndex(const vector<string>& chrom_files,
                             Genome* genome) {
  cerr << "[READING CHROMS] " << endl;
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
  for (uint32_t i = 0; i < num_of_chroms; ++i) {
    cerr << "[Chromosome " << i + 1 << "/" <<  num_of_chroms << "]" << endl;
    Chromosome& chrom = (*genome)[2 * i];
    Chromosome& chrom_rc = (*genome)[2 * i];

    chrom.name = chrom_names[i];
    chrom.length = chrom_seqs[i].size();
    chrom.strand = '+';
    chrom.sequence.resize(chrom.length);
    for (uint32_t j = 0; j < chrom.length; ++j) {
      chrom.sequence[j] = chrom_seqs[i][j];
    }
    cerr << "length = " << chrom.length << endl;
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
    BuildHashTable(&chrom);
    BuildHashTable(&chrom_rc);
  }
}

void WriteIndex(const string& index_file, const Genome& genome) {
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

    /* write hash table to disk */
    uint32_t num_of_keys = genome[i].hash_table.size();
    fwrite(&(num_of_keys), sizeof(uint32_t), 1, fout);
    for (HashTable::const_iterator it = genome[i].hash_table.begin();
        it != genome[i].hash_table.end(); ++it) {
      uint32_t hash_key = it->first;
      fwrite(&(hash_key), sizeof(uint32_t), 1, fout);
      uint32_t num_of_values = it->second.size();
      fwrite(&(num_of_values), sizeof(uint32_t), 1, fout);
      fwrite(&(it->second[0]), sizeof(uint32_t), num_of_values, fout);
    }
  }

  fclose(fout);
}

void ReadIndex(Genome* genome, const string& index_file) {
  cerr << "[READING INDEX FROM " << index_file << "]" << endl;
  FILE * fin = fopen(index_file.c_str(), "rb");
  FILE_OPEN_CHECK(fin);

  uint32_t num_of_chroms;
  FREAD_CHECK(fread(&num_of_chroms, sizeof(uint32_t), 1, fin), 1);
  genome->resize(num_of_chroms);

  char chrom_strand;
  char chrom_name[256];
  uint32_t chrom_name_len, chrom_length;
  for (uint16_t i = 0; i < num_of_chroms; ++i) {
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

    /* read hash table from disk */
    uint32_t num_of_keys = 0, hash_key = 0, num_of_values = 0;
    FREAD_CHECK(fread(&num_of_keys, sizeof(uint32_t), 1, fin), 1);
    for (uint32_t i = 0; i < num_of_keys; ++i) {
      FREAD_CHECK(fread(&hash_key, sizeof(uint32_t), 1, fin), 1);
      FREAD_CHECK(fread(&num_of_values, sizeof(uint32_t), 1, fin), 1);
      vector<uint32_t> hash_values(num_of_values);
      FREAD_CHECK(
          fread(&(hash_values[0]), sizeof(uint32_t), num_of_values, fin),
          num_of_values);

      (*genome)[i].hash_table.insert(make_pair(hash_key, hash_values));
    }
  }

  fclose(fin);
}
