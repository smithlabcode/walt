/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SMITHLAB_UTILS_HPP
#define SMITHLAB_UTILS_HPP

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <ostream>
#include <limits>
#include <iterator>
#include <cassert>
#include <cmath>


namespace smithlab {
  template <class In, class Out, class Pred> Out
  copy_if(In first, In last, Out res, Pred p) {
    while (first != last) {
      if (p(*first)) *res++ = *first;
      ++first;
    }
    return res;
  }
};

typedef size_t MASK_t;

namespace smithlab {

  // Code dealing with false discovery rate
  double get_fdr_cutoff(const size_t n_tests, std::vector<double> &pvals, 
                        const double alpha);

  void correct_pvals(const size_t n_tests, std::vector<double> &pvals); 

  // Assumes 4 nucleotide DNA alphabet
  static const size_t alphabet_size = 4;

  std::vector<std::string> split(std::string, const char *, 
				 bool get_empty_fields = false);
  std::vector<std::string> split_whitespace_quoted(std::string to_split);
  void split_whitespace(const std::string& s, std::vector<std::string> &v);
  
  std::string strip(const std::string& s);

  std::vector<std::string> 
  squash(const std::vector<std::string> &v);
  
  template <class T> std::string toa(T t) {
    std::ostringstream s;
    s << t;
    return s.str();
  }
  double log_sum_log_vec(const std::vector<double> &vals, size_t lim);
}

namespace smithlab_bits {
  static const MASK_t low_bit = 1ul;
  static const MASK_t high_bit = static_cast<size_t>(0x8000000000000000ul);
  static const MASK_t all_ones = static_cast<size_t>(-1);
  static const MASK_t all_zeros = 0ul;
  static const size_t word_size = 64ul;
}

inline size_t
base2int_upper_only(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 2;
  case 'T' : return 3;
  default  : return 4;
  }
}

inline size_t
base2int(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 2;
  case 'T' : return 3;
  case 'a' : return 0;
  case 'c' : return 1;
  case 'g' : return 2;
  case 't' : return 3; 
  default  : return 4;
  }
}

inline size_t
base2int_bs_upper_only(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 3;
  case 'G' : return 2;
  case 'T' : return 3;
  default  : return 4;
  }
}

inline size_t
base2int_bs(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 3;
  case 'G' : return 2;
  case 'T' : return 3;
  case 'a' : return 0;
  case 'c' : return 3;
  case 'g' : return 2;
  case 't' : return 3; 
  default  : return 4;
  }
}


inline size_t
base2int_bs_ag_upper_only(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 0;
  case 'T' : return 3;
  default  : return 4;
  }
}



inline size_t
base2int_bs_ag(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 0;
  case 'T' : return 3;
  case 'a' : return 0;
  case 'c' : return 1;
  case 'g' : return 0;
  case 't' : return 3; 
  default  : return 4;
  }
}


inline size_t
base2int_bs_rc(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 0;
  case 'T' : return 3;
  case 'a' : return 0;
  case 'c' : return 1;
  case 'g' : return 0;
  case 't' : return 3; 
  }
  return 4;
}

inline size_t
base2int_rc(char c) {
  switch(c) {
  case 'A' : return 3;
  case 'C' : return 2;
  case 'G' : return 1;
  case 'T' : return 0;
  case 'a' : return 3;
  case 'c' : return 2;
  case 'g' : return 1;
  case 't' : return 0; 
  }
  return 4;
}

inline char
int2base(size_t c) {
  switch(c) {
  case 0 : return 'A';
  case 1 : return 'C';
  case 2 : return 'G';
  case 3 : return 'T';
  }
  return 'N';
}

inline char
int2base_rc(size_t c) {
  switch(c) {
  case 3 : return 'A';
  case 2 : return 'C';
  case 1 : return 'G';
  case 0 : return 'T';
  }
  return 'N';
}

inline bool
isvalid(char c) {
  return (base2int(c) != 4);
}

inline std::string
i2mer(size_t n, size_t index) {
  std::string s(n, ' ');
  do {
    --n;
    s[n] = int2base(index % smithlab::alphabet_size);
    index /= smithlab::alphabet_size;
  } while (n > 0);
  return s;
}

inline std::string
i2mer_rc(size_t n, size_t index) {
  std::string s(n, ' ');
  do {
    --n;
    s[n] = int2base_rc(index % smithlab::alphabet_size);
    index /= smithlab::alphabet_size;
  } while (n > 0);
  return s;
}

inline size_t
mer2i(const std::string::const_iterator a, std::string::const_iterator b) {
  size_t multiplier = 1, index = 0;
  do {
    --b;
    index += base2int(*b)*multiplier;
    multiplier *= smithlab::alphabet_size;
  } while (b > a);
  return index;
}

inline size_t
mer2i_rc(std::string::const_iterator a, const std::string::const_iterator b) {
  size_t multiplier = 1, index = 0;
  do {
    index += base2int_rc(*a)*multiplier;
    multiplier *= smithlab::alphabet_size;
  } while (++a < b);
  return index;
}

struct SMITHLABException {
  SMITHLABException(std::string m) : message(m) {}
  std::string what() const {return message;}
  std::string message;
};


template <class T> std::string toa(T t) {
  std::ostringstream s;
  s << t;
  return s.str();
}


////////////////////////////////////////////////////////////////////////
// Code for dealing with the DNA alphabet

char
complement(int i);

inline std::string
revcomp(const std::string& s) {
  std::string r;
  std::transform(s.begin(), s.end(), back_inserter(r), complement);
  std::reverse(r.begin(), r.end());
  return r;
}

inline void
revcomp_inplace(std::string& s) {
  std::transform(s.begin(), s.end(), s.begin(), complement);
  std::reverse(s.begin(), s.end());
}

inline void
revcomp_inplace(std::string::iterator first, std::string::iterator last) {
  std::transform(first, last, first, complement);
  std::reverse(first, last);
}

inline std::string 
bits2string_masked(size_t mask, size_t bits) {
  std::string s;
  size_t selector = smithlab_bits::high_bit;
  for (size_t i = 0; i < smithlab_bits::word_size; ++i) {
    s += (selector & bits & mask) ? '1' : '0';
    selector >>= 1;
  }
  return s;
}

inline std::string 
bits2string_for_positions(size_t positions, size_t bits) {
  std::string s;
  size_t selector = smithlab_bits::high_bit;
  for (size_t i = 0; i < smithlab_bits::word_size; ++i) {
    s += (selector & bits) ? '1' : '0';
    selector >>= 1;
  }
  return s.substr(s.length() - positions);
}

inline size_t
percent(const size_t a, const size_t b) {
  return static_cast<size_t>((100.0*a)/b);
}

/////////////////////////////////////////////////////////////////////////
// We are using a conservative approach to clip that adaptors in which
// the adaptor sequence is only required to match some at some initial
// portion, and then the rest of the read is not examined.

const size_t head_length = 14;
const size_t sufficient_head_match = 11;
const size_t min_overlap = 5;

inline size_t
similarity(const std::string &s, const size_t pos, const std::string &adaptor) {
  const size_t lim = std::min(std::min(s.length() - pos, adaptor.length()), head_length);
  size_t count = 0;
  for (size_t i = 0; i < lim; ++i)
    count += (s[pos + i] == adaptor[i]);
  return count;
}

inline size_t
clip_adaptor_from_read(const std::string &adaptor, std::string &s) {
  size_t lim1 = s.length() - head_length + 1;
  for (size_t i = 0; i < lim1; ++i)
    if (similarity(s, i, adaptor) >= sufficient_head_match) {
      fill(s.begin() + i, s.end(), 'N');
      return s.length() - i;
    }
  const size_t lim2 = s.length() - min_overlap + 1;
  for (size_t i = lim1; i < lim2; ++i)
    if (similarity(s, i, adaptor) >= s.length() - i - 1) {
      fill(s.begin() + i, s.end(), 'N');
      return s.length() - i;
    }
  return 0;
}

inline void
extract_adaptors(const std::string &adaptor, std::string & T_adaptor, std::string & A_adaptor) {
  const size_t sep_idx = adaptor.find_first_of(":");
  if (adaptor.find_last_of(":") != sep_idx)
    throw SMITHLABException("ERROR: adaptor format \"T_adaptor[:A_adaptor]\"");
  if (sep_idx == std::string::npos)
    T_adaptor = A_adaptor = adaptor;
  else {
    T_adaptor = adaptor.substr(0, sep_idx);
    A_adaptor = adaptor.substr(sep_idx + 1);
  }
}
/////////////////////////////////////////////////////////////////////////

#endif
