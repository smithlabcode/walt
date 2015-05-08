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

#include <limits>
#include <string>
#include <iostream>

using std::string;
using std::ostream;
using std::cerr;
using std::endl;

#ifdef SEEDPATTERN7
#define SEEPATTERNLEN 7
#define HASHLEN 21
const uint32_t F2SEEDWIGTH = 12;
const uint32_t F2SEEDPOSITION_SIZE = 60;
const uint32_t F2SEEDPATTERN[] = {
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
    /* 15 */1, 1, 1, 0, 1, 0, 0   //
    };

const uint32_t F2SEEDPOSITION[] = {
     /* 1 */0,   1,  2, 4,  //
     /* 2 */7,   8,  9, 11,  //
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
    /* 15 */98, 99, 100, 102 //
    };
#endif

#define SEEDPATTERN3
#ifdef SEEDPATTERN3
#define SEEPATTERNLEN 3
#define HASHLEN 39
const uint32_t F2SEEDWIGTH = 13;
const uint32_t F2SEEDPATTERN[] = {
     /* 1 */0, 1, 0,  //
     /* 2 */0, 1, 0,  //
     /* 3 */0, 1, 0,  //
     /* 4 */0, 1, 0,  //
     /* 5 */0, 1, 0,  //
     /* 6 */0, 1, 0,  //
     /* 7 */0, 1, 0,  //
     /* 8 */0, 1, 0,  //
     /* 9 */0, 1, 0,  //
    /* 10 */0, 1, 0,  //
    /* 11 */0, 1, 0,  //
    /* 12 */0, 1, 0,  //
    /* 13 */0, 1, 0,  //
    /* 14 */0, 1, 0,  //
    /* 15 */0, 1, 0,  //
    /* 16 */0, 1, 0,  //
    /* 17 */0, 1, 0,  //
    /* 18 */0, 1, 0,  //
    /* 19 */0, 1, 0,  //
    /* 20 */0, 1, 0,  //
    /* 21 */0, 1, 0,  //
    /* 22 */0, 1, 0,  //
    /* 23 */0, 1, 0,  //
    /* 24 */0, 1, 0,  //
    /* 25 */0, 1, 0,  //
    /* 26 */0, 1, 0,  //
    /* 27 */0, 1, 0,  //
    /* 28 */0, 1, 0,  //
    /* 29 */0, 1, 0,  //
    /* 30 */0, 1, 0,  //
    /* 31 */0, 1, 0,  //
    /* 32 */0, 1, 0,  //
    /* 33 */0, 1, 0,  //
    /* 34 */0, 1, 0,  //
    /* 35 */0, 1, 0,  //
    /* 36 */0, 1, 0,  //
    /* 37 */0, 1, 0,  //
    /* 38 */0, 1, 0,  //
    /* 39 */0, 1, 0,  //
    /* 40 */0, 1, 0,  //
    /* 41 */0, 1, 0,  //
    /* 42 */0, 1, 0,  //
    /* 43 */0, 1, 0,  //
    /* 44 */0, 1, 0,  //
    /* 45 */0, 1, 0   //
    };

const uint32_t F2SEEDPOSITION_SIZE = 45;
const uint32_t F2SEEDPOSITION[] = {  1,   4,   7,  10,  13,  16,  19,  22,  25,  28,
                                    31,  34,  37,  40,  43,  46,  49,  52,  55,  58,
                                    61,  64,  67,  70,  73,  76,  79,  82,  85,  88,
                                    91,  94,  97, 100, 103, 106, 109, 112, 115, 118,
                                   121, 124, 127, 130, 133 };
#endif

const uint32_t MAX_UINT32 = std::numeric_limits<uint32_t>::max();

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

#ifdef DEBUG
#define DEBUG_INFO(msg, delim) { \
  cerr << msg << delim; \
}
#else
#define DEBUG_INFO(msg, delim)
#endif

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
inline uint32_t getBits(const char& nt) {
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

/* return x^p */
inline uint32_t power(const uint32_t& x, const uint32_t& p) {
  uint32_t ret = x;
  for (uint32_t i = 1; i < p; ++i) {
    ret *= x;
  }
  return ret;
}

/* transfer a k-mer to a integer number and use it as a key in the hash table */
inline uint32_t getHashValue(const char* nucleotides) {
  uint32_t hash_value = 0;
  for (uint32_t i = 0; i < F2SEEDWIGTH; ++i) {
    hash_value <<= 2;
    hash_value += getBits(nucleotides[F2SEEDPOSITION[i]]);
  }
  return hash_value;
}

#endif /* UTIL_H_ */
