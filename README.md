## WALT v1.0 ##

***WALT*** (Wildcard ALignment Tool) is a read mapping program for bisulfite sequencing in DNA methylation studies.


### Installation ###
(1) Download the source code from Github

    git clone git@github.com:smithlabcode/walt.git

(2) Build and Install
    
    cd walt
    make all
    make install


### Indexing Genome ###
    
    makedb -c <genome folder or file> -o <index file>

### Bisulfite Mapping ###

single-end reads

    walt -i <index file> -r <read files> -o <output file> [options]

paired-end reads

    walt -i <index file> -1 <mate_1 read files> -2 <mate_2 read files> -o <output file> [options]


### Mapping Options ###


| Option | Long Tag | Type | Default | Brief Description |
| :-------------: |:-------------:|:-----:|:-----:| :-----|
| -i      | -index | String | NULL |index file created by ***makedb*** command ( .dbindex) |
| -r      | -reads | String | NULL | list of single-end read files (.fastq or .fq) |
| -1      | -reads1 | String | NULL | list of paired-end read _1 files (.fastq or .fq) |
| -2      | -reads2 | String | NULL | list of paired-end read _2 files (.fastq or .fq) |
| -o      | -output | String | NULL | output file name (.sam or .mr) |
| -m      | -mismatch | Integer | 6 | maximum allowed mismatches |
| -N      | -number | Integer | 5000000 | number of reads to map at one loop |
| -a      | -ambiguous | Boolean | false | randomly output one mapped position for ambiguous reads |
| -u      | -unmapped | Boolean | false | output unmapped reads |
| -C      | -clip | String | empty | clip the specified adaptor |
| -A      | -ag-wild | Boolean | false | map using A/G bisulfite wildcards |
| -k      | -topk | Integer | 50 | maximum allowed mappings for a read in paired-end mapping |
| -L      | -fraglen | Integer | 1000 | max fragment length in paired-end mapping |
| -t      | -thread | Integer | 1 | number of threads for mapping |

To see the list of options, use "-?" or "-help".

### Examples ###

(1) **Indexing Genome**

For example, to make an index for UCSC hg19

	makedb -c hg19/ -o hg19.dbindex
   
or to make an index for chromosome 2

	makedb -c chr2.fa -o chr2.dbindex

The suffix of the index file should be '.dbindex'.
    
(2) **Bisulfite Mapping**

For example, to mapping reads to human genome hg19

	walt -i hg19.dbindex -r read_1.fq -o reads_1_mapping.sam
    
If mapping the reads from the *_2 reads file, the -A option should be set. This means that all Gs in the reads and genome are converted to As. If -A option is not set, all Cs in the reads and genome are converted to Ts.

    walt -i hg19.dbindex -r read_2.fq -A -o reads_2_mapping.sam
    
Additionally, WALT supports comma-separated list of read files. WALT produces one mapping output file for each read file. For single-end mapping, the output file names will be appended "_s1", "_s2", and so on. Notice: except the first file path, all other file paths cannot use '~'. For example, "-r ~/read_file1.fq,~/read_file2.fq" is not allowed. It should be "-r ~/read_file1.fq,/home/read_file2.fq", since linux system doesn't know it is a path except the first one.
	 
	 walt -i hg19.dbindex -r read_file1.fq,read_file2.fq,read_file3.fq -A -o reads_2_mapping.sam
	 
For paired-end reads, -1 and -2 options are used for the mate read files.
    
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -N 5000000 -o paired_reads_mapping.sam
    
Similarly, WALT supports comma-separated list of read files for paired-end mapping. WALT produces one mapping output file for each read file pair. The output file names will be appended "_p1", "_p2", and so on. One other thing to note is mate 1 and mate 2 paired files should be in the same order.

	 walt -i hg19.dbindex -1 read_file1_1.fq,read_file2_1.fq,read_file3_1.fq \ 
	                      -2 read_file1_2.fq,read_file2_2.fq,read_file3_2.fq \
	                      -o paired_reads_mapping.sam

The default number of maximum allowed mismatches is 6. The maximum allowed mismatches can be set using -m option.

    walt -i hg19.dbindex -r read_1.fq -m 4 -o reads_1_mapping.sam
    
The option -N sets the number of reads to mapping in each loop. If N is larger, the program takes large memory, especially for paired-end read mapping. If N is 1000000, both single-end and paired-end mapping take about 15 Gb memory. If N is 5000000, single-end mapping takes about 16 Gb memory, and paired-end mapping takes about 25 Gb memory. If N is set to be larger than 5000000, the program will set N to be 5000000 since when N is too large the program will take large memory but it will not be faster. The estimate memory for single-end mapping is 14 + N * (2 * rl + rnl + 16) / (1024^3) Gb, and for paired-end mapping is 14 + N * (4 * rl + 2 * rnl + 24 * k + 16) / (1024^3) Gb. N is the number of reads to map at one loop. rl is the length of reads (WALT supports mix of different lengths reads, so here rl is estimate of the average length). rnl is the length of read names. k is maximum allowed mappings for a read in paired-end mapping
    
    walt -i hg19.dbindex -r read_1.fq -N 1000000 -o reads_1_mapping.sam
    
To output the ambiguous mapped reads or unmapped reads, -u and -a options should be set. If the output format is MR, the ambigous mapped reads and unmapped reads will output to a separated file. For paired-end mapping, there will be a unmapped file and an ambiguous file for each of the mate 1 and mate 2 reads file. If the output format is SAM, uniquely mapped, ambiguous mapped, unmapped reads are output to the same file with different FLAG.
    
    walt -i hg19.dbindex -r read_1.fq -u -a -o reads_1_mapping.sam
    
To trim 3' end adaptor sequence, -C option should be set. For paired-end read mapping, using T_adaptor[:A_adaptor] format to set the adaptor. If only one adaptor seqeunce is given, the same adaptor sequence will be used for both T-rech and A-rich reads.

    walt -i hg19.dbindex -r read_1.fq -C AGATCGGAAGAGC -o reads_1_mapping.sam
    walt -i hg19.dbindex -r read_2.fq -A -C AGATCGGAAGAGC -o reads_1_mapping.sam
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -C AGATCGGAAGAGC -o paired_reads_mapping.sam
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -C AGATCGGAAGAGC:AGATCGG -o paired_mapping.sam
    
By default, WALT produces SAM format files. If you plan to get MR file, please make sure the suffix of the output file is ".mr".

	 walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -o paired_reads_mapping.sam
	 
	 
### Output Format ###

WALT supports Tab-delimited SAM and MR output formats. By default, WALT produces SAM format output files. To get MR format, the suffix of the output file should be ".mr".

**SAM Format**

* QNAME (read name)
* FLAG (bitwise FLAG)
* RNAME (chromosome name)
* POS (start position, 1-based, 0 means unmapped)
* MAPQ (255 in WALT)
* CIGAR (CIGAR string)
* RNEXT (chromosome name of the mate read)
* PNEXT (start position of the mate read)
* TLEN (observed segment length)
* SEQ (read sequence)
* QUAL (quality sequence read from fastq file)
* NM-tag (number of mismatches)

By default, WALT only outputs uniquely mapped reads. If -u or -a option is set, the ambiguous mapped or unmapped reads will be in the same SAM file with different FLAG. The FLAG 0x100 will be set for ambiguous mapped reads, and the FLAG 0x4 will be set for unmapped reads. 

**MR Format**

* RNAME (chromosome name)
* SPOS (start position, 0-based)
* EPOS (end position, 0-based)
* QNAME (read name)
* MISMATCH (number of mismatches)
* STRAND (forward or reverse strand)
* SEQ
* QUAL

If paired-end reads mapped in proper pair, the QNAME is added "FRAG:" in the beginning of the read name, the STRAND is the strand of the first mate mapped and SEQ and QUAL is merged according to their mapping positions. The overlap segment of SEQ and QUAL is from mate 1 or mate 2 and it is the one with less number of 'N' in the read sequence. MISMATCH is the sum of mismatches in mate 1 and mismatches in mate 2. If paired-end reads don't mapped in proper pair, they are treated as single-end reads.

    
### Contacts ###

***Haifeng Chen***   *haifengc@usc.edu*

***Andrew Smith***   *andrewds@usc.edu*

***Ting Chen***   *tingchen@usc.edu*
