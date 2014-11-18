/* The detail description of each function please refer to head file */

#include <math.h>

#include "reference.hpp"
#include "smithlab_os.hpp"

void ReadGenome::IdentifyChromosomes(const string& chrom_file) {
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

void ReadGenome::ReadChromosomes() {
  cerr << "[READING CHROMS] " << endl;
  vector < string > chroms;
  genome->all_chroms_len = 0;
  for (uint32_t i = 0; i < chrom_files.size(); ++i) {
    vector < string > tmp_chrom_names;
    vector < string > tmp_chroms;

    /* read chromosome name and seqeunce from the chromosome file */
    read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, tmp_chroms);

    for (uint16_t j = 0; j < tmp_chroms.size(); ++j) {
      genome->chrom_names.push_back(tmp_chrom_names[j]);
      chroms.push_back(tmp_chroms[j]);
      genome->all_chroms_len += tmp_chroms[j].size();
    }
  }

  /* copy chroms sequences to genome */
  genome->num_of_chroms = chroms.size();
  genome->chrom_seqs.resize(genome->all_chroms_len);
  cerr << "[THERE ARE " << genome->num_of_chroms << " CHROMOSOMES]" << endl;
  cerr << "[THE TAOTAL LENGTH OF ALL CHROMOSOMES IS " << genome->all_chroms_len
      << "]" << endl;
  uint32_t k = 0;
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    genome->chrom_sizes.push_back(chroms[i].size());
    genome->chrom_start_pos.push_back(k);
    for (uint32_t j = 0; j < chroms[i].size(); ++j) {
      genome->chrom_seqs[k++] = chroms[i][j];
    }
  }

  assert(k == genome->all_chroms_len);

  ToUpper();
  N2ACGT();
  C2T();
}

void ReadGenome::ToUpper() {
  cerr << "[LOWER NUCLEOTIDE TO UPPER acgt TO ACGT] " << endl;
  for (uint32_t i = 0; i < genome->all_chroms_len; ++i) {
    genome->chrom_seqs[i] = toupper(genome->chrom_seqs[i]);
  }
}

void ReadGenome::N2ACGT() {
  cerr << "[NUCLEOTIDE N TO ACGT] " << endl;
  srand (time(NULL));int
  r = 0;
  for (uint32_t i = 0; i < genome->all_chroms_len; ++i) {
    if ('N' == genome->chrom_seqs[i]) {
      r = rand() % 4;
      genome->chrom_seqs[i] = getNT(r);
    }
  }
}

void ReadGenome::C2T() {
  cerr << "[NUCLEOTIDE C TO T] " << endl;
  for (uint32_t i = 0; i < genome->all_chroms_len; ++i) {
    if ('C' == genome->chrom_seqs[i]) {
      genome->chrom_seqs[i] = 'T';
    }
  }
}
