/*
 *    This file contains several basic functions for global use.
 *
 *    Copyright (C) 2015 University of Southern California
 *                       Andrew D. Smith and Ting Chen
 *
 *    Authors: Authors: Haifeng Chen, Andrew D. Smith and Ting Chen
 *
 *    This file is part of WALT.
 *
 *    WALT is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    WALT is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with WALT.  If not, see <http://www.gnu.org/licenses/>.
 */

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

const char walt_version[] = "1.0";


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
const uint32_t NOCARED[3][100] = {
                            {  0,   2,   3,   5,   6,   8,   9,
                              11,  12,  14,  15,  17,  18,
                              20,  21,  23,  24,  26,  27,  29,  30,  32,  33,  35,  36,  38,  39,  41,  42,  44,  45,  47,  48,
                              50,  51,  53,  54,  56,  57,  59,  60,  62,  63,  65,  66,  68,  69,  71,  72,  74,  75,  77,  78,
                              80,  81,  83,  84,  86,  87,  89,  90,  92,  93,  95,  96,  98,  99, 101, 102, 104, 105, 107, 108,
                             110, 111, 113, 114, 116, 117, 119, 120, 122, 123, 125, 126, 128, 129, 131, 132, 134, 135, 137, 138 },
                       {  0,   1,   3,   4,   6,   7,   9,  10,
                              12,  13,  15,  16,  18,  19,
                              21,  22,  24,  25,  27,  28,  30,  31,  33,  34,  36,  37,  39,  40,  42,  43,  45,  46,  48,  49,
                              51,  52,  54,  55,  57,  58,  60,  61,  63,  64,  66,  67,  69,  70,  72,  73,  75,  76,  78,  79,
                              81,  82,  84,  85,  87,  88,  90,  91,  93,  94,  96,  97,  99, 100, 102, 103, 105, 106, 108, 109,
                             111, 112, 114, 115, 117, 118, 120, 121, 123, 124, 126, 127, 129, 130, 132, 133, 135, 136, 138, 139 },
                    {  0,  1,  2,   4,   5,   7,   8,  10,  11,
                              13,  14,  16,  17,  19,  20,
                              22,  23,  25,  26,  28,  29,  31,  32,  34,  35,  37,  38,  40,  41,  43,  44,  46,  47,  49,  50,
                              52,  53,  55,  56,  58,  59,  61,  62,  64,  65,  67,  68,  60,  71,  73,  74,  76,  77,  79,  80,
                              82,  83,  85,  86,  88,  89,  91,  92,  94,  95,  97,  98, 100, 101, 103, 104, 106, 107, 109, 110,
                             112, 113, 115, 116, 118, 119, 121, 122, 124, 125, 127, 128, 130, 131, 133, 134, 136, 137, 139, 140 }};

const uint32_t MAX_LINE_LENGTH = 1000;
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

//#define DEBUG
#ifdef DEBUG
#define DEBUG_INFO(msg) { \
  fprintf(stderr, "%s\n", msg); \
}
#else
#define DEBUG_INFO(msg)
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

/* return the maximal seed length for a particular read length */
inline uint32_t getSeedLength(const uint32_t& read_len) {
  return (read_len - SEEPATTERNLEN + 1) / SEEPATTERNLEN;
}

#endif /* UTIL_H_ */
