/*
 * This file contains several basic functions for global use.
 * */

#ifndef UTIL_H_
#define UTIL_H_

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <string>
#include <iostream>

using std::string;
using std::ostream;
using std::cerr;
using std::endl;

#define HASHLEN 13
#define MAX_LINE_LEN 100000

const double GB = 1024 * 1024 * 1024;

inline void MemoryAllocateCheck(void * pointer, const char * file, int line) {
  if (pointer == NULL) {
    printf("Memory allocate error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

inline void FileOpenCheck(FILE * pfile, const char * file, int line) {
  if (pfile == NULL) {
    printf("File open error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define MEMORY_ALLOCATE_CHECK(pointer)  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))

#define FREAD_CHECK(func, size) { \
  uint32_t s = func; \
  if(s != size) { \
    printf("read file error. --- %s:%s:%d\n", __FILE__, __func__, __LINE__); \
    exit(EXIT_FAILURE); \
  } \
}

#define TIME_INFO(func, msg) { \
  clock_t start_t, end_t; \
  start_t = clock(); \
  func; \
  end_t = clock(); \
  printf("[%s takes %.3lf seconds]\n", msg, \
         (double) ((end_t - start_t) / CLOCKS_PER_SEC )); \
}

/* transfer integer number to nucleotide */
inline char getNT(const int & nt) {
  switch (nt) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
  }
  return 'A';
}

/* transfer nucleotide to integer number */
inline int getBits(const char & nt) {
  switch (nt) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
  }
  return 4;
}

/* get the compliment strand nucleotide */
inline char complimentBase(char nt) {
  switch (nt) {
    case 'a':
      return ('t');
    case 'c':
      return ('g');
    case 'g':
      return ('c');
    case 't':
      return ('a');
    case 'A':
      return ('T');
    case 'C':
      return ('G');
    case 'G':
      return ('C');
    case 'T':
      return ('A');
    default:
      return ('N');
  }
}

string ReverseComplimentStrand(const string& dna_sequence) {
  string reverse_complement_sequence;
  uint32_t sequence_len = dna_sequence.size();
  for (uint32_t i = 0; i < sequence_len; ++i) {
    reverse_complement_sequence += complimentBase(dna_sequence[sequence_len - i - 1]);
  }
  return reverse_complement_sequence;
}

/* transfer a k-mer to a integer number and use it as a key in the hash table */
inline uint32_t getHashValue(const char* nucleotides) {
  uint32_t hashValue = 0;
  for (uint8_t i = 0; i < HASHLEN; ++i) {
    hashValue <<= 2;
    hashValue += getBits(nucleotides[i]);
  }
  return hashValue;
}

inline bool isFastqNameLine(uint64_t line_count) {
  return ((line_count & 3) == 0);
}

inline bool isFastqSequenceLine(uint64_t line_count) {
  return ((line_count & 3) == 1);
}

#endif /* UTIL_H_ */
