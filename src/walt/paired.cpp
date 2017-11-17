/*
 *    This file contains the functions for mapping paired-end reads.
 *    The detail description of each function please refer to head file.
 *
 *    Copyright (C) 2015 University of Southern California
 *                       Andrew D. Smith and Ting Chen
 *
 *    Authors: Haifeng Chen, Andrew D. Smith and Ting Chen
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

#include "paired.hpp"

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include <algorithm>

#include <omp.h>

uint32_t MAX(const uint32_t& a, const uint32_t& b) {
  return a > b ? a : b;
}

uint32_t MIN(const uint32_t& a, const uint32_t& b) {
  return a < b ? a : b;
}

int GetSAMFLAG(const bool& paired, const bool& paired_mapped,
               const bool& unmapped, const bool& next_unmapped, const bool& rev,
               const bool& next_rev, const bool& first, const bool& last,
               const bool& secondary_align) {
  int flag = paired ? 0x1 : 0x0;
  flag += paired_mapped ? 0x2 : 0x0;
  flag += unmapped ? 0x4 : 0x0;
  flag += next_unmapped ? 0x8 : 0x0;
  flag += rev ? 0x10 : 0x0;
  flag += next_rev ? 0x20 : 0x0;
  flag += first ? 0x40 : 0x0;
  flag += last ? 0x80 : 0x0;
  flag += secondary_align ? 0x100 : 0x0;

  return flag;
}

/* find the chromosome position in the forward strand */
void ForwardChromPosition(const uint32_t& genome_pos, const char& strand,
                          const uint32_t& chr_id, const uint32_t& read_len,
                          const Genome& genome, uint32_t& s, uint32_t& e) {
  s = genome_pos - genome.start_index[chr_id];
  s = strand == '+' ? s : genome.length[chr_id] - s - read_len;
  e = s + read_len;
}

void PairEndMapping(const string& org_read, const Genome& genome,
                    const HashTable& hash_table, const char& strand,
                    const bool& AG_WILDCARD, const uint32_t& max_mismatches,
                    const uint32_t& b, TopCandidates& top_match,
                    StatSingleReads& stat_single_reads) {
  uint32_t read_len = org_read.size();
  if (read_len < MINIMALREADLEN) {
    stat_single_reads.num_of_short_reads++;
    return;
  }

  /* return the maximal seed length for a particular read length */
  uint32_t seed_pattern_repeats = (read_len - SEEPATTERNLEN + 1)
      / SEEPATTERNLEN;
  seed_pattern_repeats = seed_pattern_repeats < 50 ? seed_pattern_repeats : 50;
  uint32_t seed_len = seed_pattern_repeats * SEEPATTERNCAREDWEIGHT;

  string read;
  if (AG_WILDCARD) {
    G2A(org_read, read_len, read);
  } else {
    C2T(org_read, read_len, read);
  }

  uint32_t cur_max_mismatches = max_mismatches;
  for (uint32_t seed_i = 0; seed_i < SEEPATTERNLEN; ++seed_i) {
    /* all exact matches are covered by the first seed */
    if (!top_match.Empty() && top_match.Full() && top_match.Top().mismatch == 0
        && seed_i)
      break;

#if defined SEEDPATTERN3 || SEEDPATTERN5
    /* all matches with 1 mismatch are covered by the first two seeds */
    if (!top_match.Empty() && top_match.Full() && top_match.Top().mismatch == 1
        && seed_i >= 2)
      break;
#endif

#ifdef SEEDPATTERN7
    /* all matches with 1 mismatch are covered by the first two seeds */
    if (!top_match.Empty() && top_match.Full() && top_match.Top().mismatch == 1
        && seed_i >= 4)
      break;
#endif

    string read_seed = read.substr(seed_i);
    uint32_t hash_value = getHashValue(read_seed.c_str());
    pair<uint32_t, uint32_t> region;
    region.first = hash_table.counter[hash_value];
    region.second = hash_table.counter[hash_value + 1];

    if (region.first == region.second)
      continue;

    IndexRegion(read_seed, genome, hash_table, seed_len, region);
    if (region.second - region.first + 1 > b) {
      continue;
    }

    for (uint32_t j = region.first; j <= region.second; ++j) {
      uint32_t genome_pos = hash_table.index[j];
      uint32_t chr_id = getChromID(genome.start_index, genome_pos);
      if (genome_pos - genome.start_index[chr_id] < seed_i)
        continue;
      genome_pos = genome_pos - seed_i;
      if (genome_pos + read_len >= genome.start_index[chr_id + 1])
        continue;

      /* check the position */
      uint32_t num_of_mismatch = 0;
      uint32_t num_of_nocared = seed_pattern_repeats * SEEPATTERNNOCAREDWEIGHT
          + seed_i;
      for (uint32_t p = 0;
          p < num_of_nocared && num_of_mismatch <= cur_max_mismatches; ++p) {
        if (genome.sequence[genome_pos + F2NOCAREDPOSITION[seed_i][p]]
            != read[F2NOCAREDPOSITION[seed_i][p]]) {
          num_of_mismatch++;
        }
      }
      for (uint32_t p = seed_pattern_repeats * SEEPATTERNLEN + seed_i;
          p < read_len && num_of_mismatch <= cur_max_mismatches; ++p) {
        if (genome.sequence[genome_pos + p] != read[p]) {
          num_of_mismatch++;
        }
      }

      if (num_of_mismatch > max_mismatches) {
        continue;
      }
      top_match.Push(CandidatePosition(genome_pos, strand, num_of_mismatch));
      if (top_match.Full()) {
        cur_max_mismatches = top_match.Top().mismatch;
      }
    }
  }
}

void OutputPairedStatInfo(const StatPairedReads& stat_paired_reads,
                          const string& output_file) {
  FILE * fstat = fopen(string(output_file + ".mapstats").c_str(), "w");
  fprintf(fstat, "[TOTAL NUMBER OF READ PAIRS: %u]\n",
          stat_paired_reads.total_read_pairs);
  fprintf(
      fstat,
      "[UNIQUELY MAPPED READ PAIRS: %u (%.2lf%%)]\n",
      stat_paired_reads.unique_mapped_pairs,
      100.00 * stat_paired_reads.unique_mapped_pairs
          / stat_paired_reads.total_read_pairs);
  fprintf(
      fstat,
      "[AMBIGUOUS MAPPED READ PAIRS: %u (%.2lf%%)]\n",
      stat_paired_reads.ambiguous_mapped_pairs,
      100.00 * stat_paired_reads.ambiguous_mapped_pairs
          / stat_paired_reads.total_read_pairs);
  fprintf(
      fstat,
      "[UNMAPPED READS PAIRS: %u (%.2lf%%)]\n",
      stat_paired_reads.unmapped_pairs,
      100.00 * stat_paired_reads.unmapped_pairs
          / stat_paired_reads.total_read_pairs);
  //////////////////MATE 1////////////////////////
  fprintf(
      fstat,
      "   [UNIQUELY MAPPED READS IN MATE_1: %u (%.2lf%%)]\n",
      stat_paired_reads.stat_single_reads_1.unique_mapped_reads,
      100.00 * stat_paired_reads.stat_single_reads_1.unique_mapped_reads
          / stat_paired_reads.total_read_pairs);
  fprintf(
      fstat,
      "   [AMBIGUOUS MAPPED READS IN MATE_1: %u (%.2lf%%)]\n",
      stat_paired_reads.stat_single_reads_1.ambiguous_mapped_reads,
      100.00 * stat_paired_reads.stat_single_reads_1.ambiguous_mapped_reads
          / stat_paired_reads.total_read_pairs);
  fprintf(
      fstat,
      "   [UNMAPPED READS IN MATE_1: %u (%.2lf%%)]\n",
      stat_paired_reads.stat_single_reads_1.unmapped_reads,
      100.00 * stat_paired_reads.stat_single_reads_1.unmapped_reads
          / stat_paired_reads.total_read_pairs);
  //////////////////MATE 2////////////////////////
  fprintf(
      fstat,
      "   [UNIQUELY MAPPED READS IN MATE_2: %u (%.2lf%%)]\n",
      stat_paired_reads.stat_single_reads_2.unique_mapped_reads,
      100.00 * stat_paired_reads.stat_single_reads_2.unique_mapped_reads
          / stat_paired_reads.total_read_pairs);
  fprintf(
      fstat,
      "   [AMBIGUOUS MAPPED READS IN MATE_2: %u (%.2lf%%)]\n",
      stat_paired_reads.stat_single_reads_2.ambiguous_mapped_reads,
      100.00 * stat_paired_reads.stat_single_reads_2.ambiguous_mapped_reads
          / stat_paired_reads.total_read_pairs);
  fprintf(
      fstat,
      "   [UNMAPPED READS IN MATE_2: %u (%.2lf%%)]\n",
      stat_paired_reads.stat_single_reads_2.unmapped_reads,
      100.00 * stat_paired_reads.stat_single_reads_2.unmapped_reads
          / stat_paired_reads.total_read_pairs);

  if (stat_paired_reads.stat_single_reads_1.num_of_short_reads != 0
      || stat_paired_reads.stat_single_reads_2.num_of_short_reads != 0) {
    fprintf(fstat, "\n\n[READS SHORTER THAN %d ARE IGNORED]\n",
            MINIMALREADLEN);
    fprintf(
        fstat,
        "[%u (%.2lf%%) READS ARE SHORTER THAN %d IN MATE_1]\n",
        stat_paired_reads.stat_single_reads_1.num_of_short_reads,
        100.00 * stat_paired_reads.stat_single_reads_1.num_of_short_reads
            / stat_paired_reads.total_read_pairs,
        MINIMALREADLEN);
    fprintf(
        fstat,
        "[%u (%.2lf%%) READS ARE SHORTER THAN %d IN MATE_2]\n",
        stat_paired_reads.stat_single_reads_2.num_of_short_reads,
        100.00 * stat_paired_reads.stat_single_reads_2.num_of_short_reads
            / stat_paired_reads.total_read_pairs,
        MINIMALREADLEN);
  }

  // distribution of fragment length
  double average_fragment_len = 0.0;
  fprintf(fstat, "\n\nDISTRIBUTION OF PAIRED-END FRAGMNET LENGTH\n");
  for (uint32_t i = 1; i < stat_paired_reads.fragment_len_count.size(); ++i) {
    if (stat_paired_reads.fragment_len_count[i] != 0) {
      fprintf(
          fstat,
          "%u:\t%u\t(%.2lf%%)\n",
          i,
          stat_paired_reads.fragment_len_count[i],
          100.00 * stat_paired_reads.fragment_len_count[i]
              / stat_paired_reads.unique_mapped_pairs);
      average_fragment_len += i * stat_paired_reads.fragment_len_count[i];
    }
  }
  average_fragment_len /= stat_paired_reads.unique_mapped_pairs;

  // standard deviation
  double standard_deviation = 0.0;
  for (uint32_t i = 1; i < stat_paired_reads.fragment_len_count.size(); ++i) {
    if (stat_paired_reads.fragment_len_count[i] != 0) {
      standard_deviation += (average_fragment_len - i)
          * (average_fragment_len - i)
          * stat_paired_reads.fragment_len_count[i];
    }
  }
  standard_deviation = sqrt(standard_deviation);
  fprintf(fstat, "\n\nAVERAGE VALUE OF PAIRED-END FRAGMNET LENGTHS: %lf\n",
          average_fragment_len);
  fprintf(fstat, "STANDARD DEVIATION OF PAIRED-END FRAGMNET LENGTHS: %lf\n",
          standard_deviation);

  fclose(fstat);
}

int OutputBestPairedResults(const CandidatePosition& r1,
                            const CandidatePosition& r2, const int& frag_range,
                            const uint32_t& read_len1,
                            const uint32_t& read_len2, const Genome& genome,
                            const string& read_name, const string& read_seq1,
                            const string& read_score1, const string& read_seq2,
                            const string& read_score2, const bool& SAM,
                            const bool& PBAT,
                            FILE * fout) {
  string read_seq2_rev = ReverseComplimentString(read_seq2);
  string read_scr2_rev = ReverseString(read_score2);

  uint32_t chr_id1 = getChromID(genome.start_index, r1.genome_pos);
  uint32_t chr_id2 = getChromID(genome.start_index, r2.genome_pos);

  uint32_t s1 = 0, s2 = 0, e1 = 0, e2 = 0;
  ForwardChromPosition(r1.genome_pos, r1.strand, chr_id1, read_len1, genome, s1,
                       e1);
  ForwardChromPosition(r2.genome_pos, r2.strand, chr_id2, read_len2, genome, s2,
                       e2);

  uint32_t overlap_s = MAX(s1, s2);
  uint32_t overlap_e = MIN(e1, e2);

  uint32_t one_l = r1.strand == '+' ? s1 : MAX(overlap_e, s1);
  uint32_t one_r = r1.strand == '+' ? MIN(overlap_s, e1) : e1;

  uint32_t two_l = r1.strand == '+' ? MAX(overlap_e, s2) : s2;
  uint32_t two_r = r1.strand == '+' ? e2 : MIN(overlap_s, e2);

  int len = r1.strand == '+' ? (two_r - one_l) : (one_r - two_l);
  if (SAM) {
    return len;
  }

  string seq(len, 'N');
  string scr(len, 'B');
  if (len > 0 && len <= frag_range) {
    // lim_one: offset in merged sequence where overlap starts
    uint32_t lim_one = one_r - one_l;
    copy(read_seq1.begin(), read_seq1.begin() + lim_one, seq.begin());
    copy(read_score1.begin(), read_score1.begin() + lim_one, scr.begin());

    uint32_t lim_two = two_r - two_l;
    copy(read_seq2_rev.end() - lim_two, read_seq2_rev.end(),
         seq.end() - lim_two);
    copy(read_scr2_rev.end() - lim_two, read_scr2_rev.end(),
         scr.end() - lim_two);

    // deal with overlapping part
    if (overlap_s < overlap_e) {
      uint32_t one_bads = count(read_seq1.begin(), read_seq1.end(), 'N');
      int info_one = read_len1 - (one_bads + r1.mismatch);

      uint32_t two_bads = count(read_seq2_rev.begin(), read_seq2_rev.end(),
                                'N');
      int info_two = read_len2 - (two_bads + r2.mismatch);

      // use the mate with the most info to fill in the overlap
      if (info_one >= info_two) {
        uint32_t a = r1.strand == '+' ? (overlap_s - s1) : (e1 - overlap_e);
        uint32_t b = r1.strand == '+' ? (overlap_e - s1) : (e1 - overlap_s);
        copy(read_seq1.begin() + a, read_seq1.begin() + b,
             seq.begin() + lim_one);
        copy(read_score1.begin() + a, read_score1.begin() + b,
             scr.begin() + lim_one);
      } else {
        uint32_t a = r1.strand == '+' ? (overlap_s - s2) : (e2 - overlap_e);
        uint32_t b = r1.strand == '+' ? (overlap_e - s2) : (e2 - overlap_s);
        copy(read_seq2_rev.begin() + a, read_seq2_rev.begin() + b,
             seq.begin() + lim_one);
        copy(read_scr2_rev.begin() + a, read_scr2_rev.begin() + b,
             scr.begin() + lim_one);
      }
    }
  }

  char strand = r1.strand;
  uint32_t start_pos = r1.strand == '+' ? s1 : s2;
  if (PBAT) {
    // if PBAT, the whole fragment is considered as A-rich;
    // therefore it is necessary to report a reverse-complement
    strand = r1.strand == '+' ? '-' : '+';
    seq = ReverseComplimentString(seq);
    scr = ReverseString(scr);
  }
  fprintf(fout, "%s\t%u\t%u\tFRAG:%s\t%u\t%c\t%s\t%s\n",
          genome.name[chr_id1].c_str(), start_pos, start_pos + len,
          read_name.c_str(), r1.mismatch + r2.mismatch, strand, seq.c_str(),
          scr.c_str());

  return len;
}

void GetBestMatch4Single(const vector<CandidatePosition>& ranked_results,
                         const int ranked_results_size, const Genome& genome,
                         const uint32_t& read_len, const string& read_name,
                         const string& read_seq, const string& read_score,
                         const uint32_t& max_mismatches,
                         BestMatch& best_match) {
  for (int i = ranked_results_size - 1; i >= 0; --i) {
    const CandidatePosition& r = ranked_results[i];
    if (r.mismatch < best_match.mismatch) {
      best_match = BestMatch(r.genome_pos, 1, r.strand, r.mismatch);
    } else if (r.mismatch == best_match.mismatch) {
      if (best_match.genome_pos == r.genome_pos) {
        continue;
      } else {
        best_match.genome_pos = r.genome_pos;
        best_match.strand = r.strand;
        best_match.times++;
      }
    } else {
      break;
    }
  }
}

int GetFragmentLength(const CandidatePosition& r1, const CandidatePosition& r2,
                      const uint32_t& frag_range, const uint32_t& read_len1,
                      const uint32_t& read_len2, const Genome& genome,
                      const uint32_t& chr_id1, const uint32_t& chr_id2) {
  uint32_t s1 = 0, s2 = 0, e1 = 0, e2 = 0;
  ForwardChromPosition(r1.genome_pos, r1.strand, chr_id1, read_len1, genome, s1,
                       e1);
  ForwardChromPosition(r2.genome_pos, r2.strand, chr_id2, read_len2, genome, s2,
                       e2);

  return r1.strand == '+' ? (e2 - s1) : (e1 - s2);
}

void OutputPairedSAM(const BestMatch& best_match_1,
                     const BestMatch& best_match_2, const Genome& genome,
                     const string& read_name, const string& read_seq1,
                     const string& read_score1, const string& read_seq2,
                     const string& read_score2, const int& len,
                     const int& flag_1, const int& flag_2,
                     StatPairedReads& stat_paired_reads, FILE * fout) {
  uint32_t chr_id_1 = getChromID(genome.start_index, best_match_1.genome_pos);
  uint32_t chr_id_2 = getChromID(genome.start_index, best_match_2.genome_pos);
  uint32_t s1 = 0, s2 = 0, e1 = 0, e2 = 0;
  ForwardChromPosition(best_match_1.genome_pos, best_match_1.strand, chr_id_1,
                       read_seq1.size(), genome, s1, e1);
  ForwardChromPosition(best_match_2.genome_pos, best_match_2.strand, chr_id_2,
                       read_seq2.size(), genome, s2, e2);

  uint32_t mismatch1 = best_match_1.mismatch;
  uint32_t mismatch2 = best_match_2.mismatch;
  if (best_match_1.times == 0) {
    s1 = 0;
    mismatch1 = 0;
  } else {
    s1 += 1;
  }
  if (best_match_2.times == 0) {
    s2 = 0;
    mismatch2 = 0;
  } else {
    s2 += 1;
  }

  int len1 = best_match_1.strand == '+' ? len : -len;
  int len2 = best_match_2.strand == '+' ? len : -len;

  string rnext1 = "=", rnext2 = "=";
  if (flag_1 & 0x2) {
    rnext1 = "=";
    rnext2 = "=";
  } else {
    if (best_match_1.times == 0) {
      rnext1 = "*";
    } else {
      rnext1 = genome.name[chr_id_1].c_str();
    }

    if (best_match_2.times == 0) {
      rnext2 = "*";
    } else {
      rnext2 = genome.name[chr_id_2].c_str();
    }
  }

  string read_seq1_tmp = read_seq1;
  string read_seq2_tmp = read_seq2;
  string read_score1_tmp = read_score1;
  string read_score2_tmp = read_score2;
  if (best_match_1.strand == '-') {
    read_seq1_tmp = ReverseComplimentString(read_seq1_tmp);
    read_score1_tmp = ReverseString(read_score1_tmp);
  }
  if (best_match_2.strand == '-') {
    read_seq2_tmp = ReverseComplimentString(read_seq2_tmp);
    read_score2_tmp = ReverseString(read_score2_tmp);
  }

  uint32_t read_len1 = read_seq1.size();
  uint32_t read_len2 = read_seq2.size();

  if (best_match_1.times == 0
      && stat_paired_reads.stat_single_reads_1.unmapped) {
    fprintf(fout, "%s\t%d\t*\t%u\t255\t*\t%s\t%u\t%d\t%s\t%s\tNM:i:%u\n",
            read_name.c_str(), flag_1, s1, rnext2.c_str(), s2, len1,
            read_seq1_tmp.c_str(), read_score1_tmp.c_str(), mismatch1);
  } else if (best_match_1.times == 1) {
    fprintf(fout, "%s\t%d\t%s\t%u\t255\t%uM\t%s\t%u\t%d\t%s\t%s\tNM:i:%u\n",
            read_name.c_str(), flag_1, genome.name[chr_id_1].c_str(), s1,
            read_len1, rnext2.c_str(), s2, len1, read_seq1_tmp.c_str(),
            read_score1_tmp.c_str(), mismatch1);
  } else if (best_match_1.times >= 2
      && stat_paired_reads.stat_single_reads_1.ambiguous) {
    fprintf(fout, "%s\t%d\t%s\t%u\t255\t%uM\t%s\t%u\t%d\t%s\t%s\tNM:i:%u\n",
            read_name.c_str(), flag_1, genome.name[chr_id_1].c_str(), s1,
            read_len1, rnext2.c_str(), s2, len1, read_seq1_tmp.c_str(),
            read_score1_tmp.c_str(), mismatch1);
  }

  if (best_match_2.times == 0
      && stat_paired_reads.stat_single_reads_2.unmapped) {
    fprintf(fout, "%s\t%d\t*\t%u\t255\t*\t%s\t%u\t%d\t%s\t%s\tNM:i:%u\n",
            read_name.c_str(), flag_2, s2, rnext1.c_str(), s1, len2,
            read_seq2_tmp.c_str(), read_score2_tmp.c_str(), mismatch2);
  } else if (best_match_2.times == 1) {
    fprintf(fout, "%s\t%d\t%s\t%u\t255\t%uM\t%s\t%u\t%d\t%s\t%s\tNM:i:%u\n",
            read_name.c_str(), flag_2, genome.name[chr_id_2].c_str(), s2,
            read_len2, rnext1.c_str(), s1, len2, read_seq2_tmp.c_str(),
            read_score2_tmp.c_str(), mismatch2);
  } else if (best_match_2.times >= 2
      && stat_paired_reads.stat_single_reads_2.ambiguous) {
    fprintf(fout, "%s\t%d\t%s\t%u\t255\t%uM\t%s\t%u\t%d\t%s\t%s\tNM:i:%u\n",
            read_name.c_str(), flag_2, genome.name[chr_id_2].c_str(), s2,
            read_len2, rnext1.c_str(), s1, len2, read_seq2_tmp.c_str(),
            read_score2_tmp.c_str(), mismatch2);
  }
}

/* merge the mapping results from paired reads */
void MergePairedEndResults(
    const Genome& genome, const string& read_name, const string& read_seq1,
    const string& read_score1, const string& read_seq2,
    const string& read_score2,
    const vector<vector<CandidatePosition> >& ranked_results,
    const vector<int>& ranked_results_size, const int& frag_range,
    const uint32_t& max_mismatches, const bool& SAM,
    const bool& PBAT,
    StatPairedReads& stat_paired_reads, FILE * fout) {
#ifdef DEBUG
  for (int i = ranked_results_size[0] - 1; i >= 0; --i) {
    const CandidatePosition& r1 = ranked_results[0][i];
    uint32_t chr_id1 = getChromID(genome.start_index, r1.genome_pos);
    uint32_t start_pos = r1.genome_pos - genome.start_index[chr_id1];
    if ('-' == r1.strand) {
      start_pos = genome.length[chr_id1] - start_pos - read_seq1.size();
    }
    uint32_t end_pos = start_pos + read_seq1.size();
    fprintf(stderr, "%u %s %u %u %c %u\n", r1.genome_pos,
            genome.name[chr_id1].c_str(), start_pos, end_pos, r1.strand,
            r1.mismatch);
  }
  for (int j = ranked_results_size[1] - 1; j >= 0; --j) {
    const CandidatePosition& r2 = ranked_results[1][j];
    uint32_t chr_id2 = getChromID(genome.start_index, r2.genome_pos);
    uint32_t start_pos = r2.genome_pos - genome.start_index[chr_id2];
    if ('-' == r2.strand) {
      start_pos = genome.length[chr_id2] - start_pos - read_seq2.size();
    }
    uint32_t end_pos = start_pos + read_seq2.size();
    fprintf(stderr, "%u %s %u %u %c %u\n", r2.genome_pos,
            genome.name[chr_id2].c_str(), start_pos, end_pos, r2.strand,
            r2.mismatch);
  }
#endif
  uint32_t read_len1 = read_seq1.size();
  uint32_t read_len2 = read_seq2.size();
  pair<int, int> best_pair(-1, -1);
  uint32_t min_num_of_mismatch = max_mismatches;
  uint64_t best_pos = 0;
  uint32_t best_times = 0;
  for (int i = ranked_results_size[0] - 1; i >= 0; --i) {
    for (int j = ranked_results_size[1] - 1; j >= 0; --j) {
      const CandidatePosition& r1 = ranked_results[0][i];
      const CandidatePosition& r2 = ranked_results[1][j];
      if (r1.strand == r2.strand)
        continue;

      uint32_t num_of_mismatch = r1.mismatch + r2.mismatch;
      if (num_of_mismatch > min_num_of_mismatch)
        break;

      uint32_t chr_id1 = getChromID(genome.start_index, r1.genome_pos);
      uint32_t chr_id2 = getChromID(genome.start_index, r2.genome_pos);
      if (chr_id1 != chr_id2)
        continue;

      int frag_size = GetFragmentLength(r1, r2, frag_range, read_len1,
                                        read_len2, genome, chr_id1, chr_id2);
      if (frag_size <= 0 || frag_size > frag_range)
        continue;

      uint64_t cur_pos = r1.genome_pos;
      cur_pos <<= 32;
      cur_pos += r2.genome_pos;
      if (num_of_mismatch < min_num_of_mismatch) {
        best_pair = make_pair(i, j);
        best_times = 1;
        min_num_of_mismatch = num_of_mismatch;
        best_pos = cur_pos;
      } else if (num_of_mismatch == min_num_of_mismatch
          && cur_pos != best_pos) {
        best_pair = make_pair(i, j);
        best_times++;
      }
    }
  }

  BestMatch best_match_1(0, 0, '+', max_mismatches);
  BestMatch best_match_2(0, 0, '+', max_mismatches);
  bool is_paired_mapped = false;
  int len = 0;
  if (best_times == 1) {
    stat_paired_reads.unique_mapped_pairs++;
    len = OutputBestPairedResults(ranked_results[0][best_pair.first],
                                  ranked_results[1][best_pair.second],
                                  frag_range, read_len1, read_len2, genome,
                                  read_name, read_seq1, read_score1, read_seq2,
                                  read_score2, SAM, PBAT, fout);
    stat_paired_reads.fragment_len_count[len]++;
    if (SAM) {  // SAM
      is_paired_mapped = true;
      const CandidatePosition& r1 = ranked_results[0][best_pair.first];
      const CandidatePosition& r2 = ranked_results[1][best_pair.second];
      best_match_1 = BestMatch(r1.genome_pos, 1, r1.strand, r1.mismatch);
      best_match_2 = BestMatch(r2.genome_pos, 1, r2.strand, r2.mismatch);
    }
  } else {
    if (best_times >= 2) {
      stat_paired_reads.ambiguous_mapped_pairs++;
    } else {
      stat_paired_reads.unmapped_pairs++;
    }
    GetBestMatch4Single(ranked_results[0], ranked_results_size[0], genome,
                        read_len1, read_name, read_seq1, read_score1,
                        max_mismatches, best_match_1);
    GetBestMatch4Single(ranked_results[1], ranked_results_size[1], genome,
                        read_len2, read_name, read_seq2, read_score2,
                        max_mismatches, best_match_2);
    StatInfoUpdate(best_match_1.times, stat_paired_reads.stat_single_reads_1);
    StatInfoUpdate(best_match_2.times, stat_paired_reads.stat_single_reads_2);
    if (!SAM) {
      OutputSingleResults(best_match_1, read_name, read_seq1, read_score1,
                          genome, PBAT, stat_paired_reads.stat_single_reads_1,
                          fout);
      OutputSingleResults(best_match_2, read_name, read_seq2, read_score2,
                          genome, !PBAT, stat_paired_reads.stat_single_reads_2,
                          fout);
    }
  }
  if (SAM) {  // Output SAM
    int flag_1 = GetSAMFLAG(true, is_paired_mapped, best_match_1.times == 0,
                            best_match_2.times == 0, best_match_1.strand == '-',
                            best_match_2.strand == '-', true, false,
                            best_match_1.times >= 2);
    int flag_2 = GetSAMFLAG(true, is_paired_mapped, best_match_2.times == 0,
                            best_match_1.times == 0, best_match_2.strand == '-',
                            best_match_1.strand == '-', false, true,
                            best_match_2.times >= 2);
    OutputPairedSAM(best_match_1, best_match_2, genome, read_name, read_seq1,
                    read_score1, read_seq2, read_score2, len, flag_1, flag_2,
                    stat_paired_reads, fout);
  }
}

void ProcessPairedEndReads(const string& command, const string& index_file,
                           const string& reads_file_p1,
                           const string& reads_file_p2,
                           const string& output_file,
                           const uint32_t& n_reads_to_process,
                           const uint32_t& max_mismatches, const uint32_t& b,
                           const string& adaptor, const bool& PBAT,
                           const uint32_t& top_k, const int& frag_range,
                           const bool& ambiguous, const bool& unmapped,
                           const bool& SAM, const int& num_of_threads) {
  // LOAD THE INDEX HEAD INFO
  Genome genome;
  HashTable hash_table;

  uint32_t size_of_index;
  ReadIndexHeadInfo(index_file, genome, size_of_index);
  genome.sequence.resize(genome.length_of_genome);
  hash_table.counter.resize(power(4, F2SEEDKEYWIGTH) + 1);
  hash_table.index.resize(size_of_index);

  vector<vector<string> > index_names(2, vector<string>(2));
  if (!PBAT) {
    index_names[0][0] = index_file + "_CT00";
    index_names[0][1] = index_file + "_CT01";
    index_names[1][0] = index_file + "_GA10";
    index_names[1][1] = index_file + "_GA11";
  } else {
    index_names[0][0] = index_file + "_GA10";
    index_names[0][1] = index_file + "_GA11";
    index_names[1][0] = index_file + "_CT00";
    index_names[1][1] = index_file + "_CT01";
  }

  vector<vector<string> > read_names(2, vector<string>(n_reads_to_process));
  vector<vector<string> > read_seqs(2, vector<string>(n_reads_to_process));
  vector<vector<string> > read_scores(2, vector<string>(n_reads_to_process));

  vector<int> ranked_results_size(2);
  vector<vector<CandidatePosition> > ranked_results(2,
          vector<CandidatePosition>(top_k));

  vector<vector<TopCandidates> > top_results(2,
         vector<TopCandidates>(n_reads_to_process));

  FILE * fin[2];
  fin[0] = fopen(reads_file_p1.c_str(), "r");
  if (!fin[0]) {
    throw SMITHLABException("cannot open input file " + reads_file_p1);
  }
  fin[1] = fopen(reads_file_p2.c_str(), "r");
  if (!fin[1]) {
    throw SMITHLABException("cannot open input file " + reads_file_p2);
  }

  string adaptors[2];
  extract_adaptors(adaptor, adaptors[0], adaptors[1]);
  clock_t start_t = clock();
  FILE * fout = fopen(output_file.c_str(), "w");
  if (!fout) {
    throw SMITHLABException("cannot open input file " + output_file);
  }
  uint32_t num_of_reads[2];
  StatPairedReads stat_paired_reads(ambiguous, unmapped, frag_range,
                                    output_file, SAM);
  bool AG_WILDCARD = true;
  fprintf(stderr, "[MAPPING PAIRED-END READS FROM THE FOLLOWING TWO FILES]\n");
  fprintf(stderr, "   %s (AND)\n   %s\n", reads_file_p1.c_str(),
          reads_file_p2.c_str());
  fprintf(stderr, "[OUTPUT MAPPING RESULTS TO %s]\n", output_file.c_str());
  if (SAM) {
    SAMHead(index_file, command, fout);
  }
  omp_set_dynamic(0);
  omp_set_num_threads(num_of_threads);
  for (uint32_t i = 0;; i += n_reads_to_process) {
    num_of_reads[0] = num_of_reads[1] = 0;
    for (uint32_t pi = 0; pi < 2; ++pi) {  // paired end reads _1 and _2
      if (!PBAT)
        AG_WILDCARD = pi == 1 ? true : false;
      else
        AG_WILDCARD = pi == 0 ? true : false;
      StatSingleReads& stat_single_reads =
          pi == 0 ?
              stat_paired_reads.stat_single_reads_1 :
              stat_paired_reads.stat_single_reads_2;
      LoadReadsFromFastqFile(fin[pi], i, n_reads_to_process, adaptors[pi],
                             num_of_reads[pi], read_names[pi], read_seqs[pi],
                             read_scores[pi]);
      if (num_of_reads[pi] == 0)
        break;

      //Initialize the paired results
      for (uint32_t j = 0; j < num_of_reads[pi]; ++j) {
        top_results[pi][j].Clear();
        top_results[pi][j].SetSize(top_k);
      }

      for (uint32_t fi = 0; fi < 2; ++fi) {
        ReadIndex(index_names[pi][fi], genome, hash_table);
        char strand = fi == 0 ? '+' : '-';

#pragma omp parallel for
        for (uint32_t j = 0; j < num_of_reads[pi]; ++j) {
          PairEndMapping(read_seqs[pi][j], genome, hash_table, strand,
                         AG_WILDCARD, max_mismatches, b, top_results[pi][j],
                         stat_single_reads);
        }
#pragma omp barrier
      }
    }
    if (num_of_reads[0] != num_of_reads[1]) {
      fprintf(stderr,
              "The number of reads in paired-end files should be the same.\n");
      exit( EXIT_FAILURE);
    }
    if (num_of_reads[0] == 0) {
      break;
    }
    stat_paired_reads.total_read_pairs += num_of_reads[0];
    ///////////////////////////////////////////////////////////
    // Merge Paired-end results
    for (uint32_t j = 0; j < num_of_reads[0]; ++j) {
      for (uint32_t pi = 0; pi < 2; ++pi) {
        ranked_results_size[pi] = 0;
        while (!top_results[pi][j].candidates.empty()) {
          ranked_results[pi][ranked_results_size[pi]++] =
              top_results[pi][j].Top();
          top_results[pi][j].Pop();
        }
      }

      MergePairedEndResults(genome, read_names[0][j], read_seqs[0][j],
                            read_scores[0][j], read_seqs[1][j],
                            read_scores[1][j], ranked_results,
                            ranked_results_size, frag_range, max_mismatches,
                            SAM, PBAT, stat_paired_reads, fout);
    }

    if (num_of_reads[0] < n_reads_to_process)
      break;
  }

  fclose(fin[0]);
  fclose(fin[1]);
  fclose(fout);

  OutputPairedStatInfo(stat_paired_reads, output_file);
  fprintf(stderr, "[MAPPING TAKES %.0lf SECONDS]\n",
          (double(clock() - start_t) / CLOCKS_PER_SEC));
}
