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

#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <limits>
#include <string>

#include "seedpattern.hpp"

const char walt_version[] = "1.0";

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
    default:
      fprintf(stderr, "[ERROR: NON-ACGT NUCLEOTIDE]");
      exit(EXIT_FAILURE);
  }
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
    default:
      fprintf(stderr, "[ERROR: NON-ACGT NUCLEOTIDE]");
      exit(EXIT_FAILURE);
  }
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
      fprintf(stderr, "[ERROR: NON-ACGT NUCLEOTIDE]");
      exit(EXIT_FAILURE);
  }
}

/* not A, C, G, or T */
inline bool nonACGT(const char& nt) {
  return !(nt == 'A' || nt == 'C' || nt == 'G' || nt == 'T');
}

/* non-ACGT nucleotide to A, C, G, or T */
inline char toACGT(const char& nt) {
  if (nonACGT(nt)) {
    int r = rand() % 4;
    return getNT(r);
  }

  return nt;
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
  for (uint32_t i = 0; i < F2SEEDKEYWIGTH; ++i) {
    hash_value <<= 2;
    hash_value += getBits(nucleotides[F2CAREDPOSITION[i]]);
  }
  return hash_value;
}

#endif /* UTIL_H_ */
