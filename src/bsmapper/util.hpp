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

#define HASHLEN 21
const uint32_t F2SEEDWIGTH = 12;

const uint32_t F2SEEDPATTERNx[] = {
/* 1 */1, 1, 1, 0, 1, 0, 0,  //
    /* 2 */1, 1, 1, 0, 1, 0, 0,  //
    /* 3 */1, 1, 1, 0, 1, 0, 0,  //
    /* 4 */1, 1, 1, 0, 1, 0, 0,  //
    /* 5 */1, 1, 1, 0, 1, 0, 0,  //
    /* 6 */1, 1, 1, 0, 1, 0, 0,  //
    /* 7 */1, 1, 1, 0, 1, 0, 0,  //
    /* 8 */1, 1, 1, 0, 1, 0, 0,  //
    /* 9 */1, 1, 1, 0, 1, 0, 0,  //
    /* 10 */1, 1, 1, 0, 1, 0, 0,  //
    /* 11 */1, 1, 1, 0, 1, 0, 0,  //
    /* 12 */1, 1, 1, 0, 1, 0, 0,  //
    /* 13 */1, 1, 1, 0, 1, 0, 0,  //
    /* 14 */1, 1, 1, 0, 1, 0, 0,  //
    };

const uint32_t F2SEEDPAOSITION[] = {
/* 1 */0, 1, 2, 4,  //
    /* 2 */7, 8, 9, 11,  //
    /* 3 */14, 15, 16, 18,  //
    /* 4 */21, 22, 23, 25,  //
    /* 5 */28, 29, 30, 32,  //
    /* 6 */35, 36, 37, 39,  //
    /* 7 */42, 43, 44, 46,  //
    /* 8 */49, 50, 51, 53,  //
    /* 9 */56, 57, 58, 60,  //
    /* 10 */63, 64, 65, 67,  //
    /* 11 */70, 71, 72, 74,  //
    /* 12 */77, 78, 79, 81,  //
    /* 13 */84, 85, 86, 88,  //
    /* 14 */91, 92, 93, 95,  //
    };

#define MAX_LINE_LEN 100000
const double GB = 1024 * 1024 * 1024;

inline void MemoryAllocateCheck(void* pointer, const char* file, int line) {
  if (pointer == NULL) {
    fprintf(stderr, "Memory allocate error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

inline void FileOpenCheck(FILE* pfile, const char* file, int line) {
  if (pfile == NULL) {
    fprintf(stderr, "File open error in %s at line %d\n", file, line);
    exit(EXIT_FAILURE);
  }
}

#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define MEMORY_ALLOCATE_CHECK(pointer)  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))

#define FREAD_CHECK(func, size) { \
  uint32_t s = func; \
  if(s != size) { \
    fprintf(stderr, "read file error. --- %s:%s:%d\n", __FILE__, __func__, __LINE__); \
    exit(EXIT_FAILURE); \
  } \
}

#define TIME_INFO(func, msg) { \
  clock_t start_t, end_t; \
  start_t = clock(); \
  func; \
  end_t = clock(); \
  fprintf(stderr, "[%s TAKES %.3lf SECONDS]\n", msg, \
         (double) ((end_t - start_t) / CLOCKS_PER_SEC )); \
}

/* transfer integer number to nucleotide */
inline char getNT(const int& nt) {
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
inline int getBits(const char& nt) {
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
inline char complimentBase(const char& nt) {
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

/* transfer a k-mer to a integer number and use it as a key in the hash table */
inline uint32_t getHashValue(const char* nucleotides) {
  uint32_t hash_value = 0;
  for (uint32_t i = 0; i < F2SEEDWIGTH; ++i) {
    hash_value <<= 2;
    hash_value += getBits(nucleotides[F2SEEDPAOSITION[i]]);
  }
  return hash_value;
}

#endif /* UTIL_H_ */
