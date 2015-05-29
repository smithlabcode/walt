## Pre-release ##

**Please download the pre-release verson for data analysis.** 


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
| -o      | -output | String | NULL | output file name |
| -m      | -mismatch | Integer | 6 | maximum allowed mismatches |
| -N      | -number | Integer | 5000000 | number of reads to map at one loop |
| -a      | -ambiguous | Boolean | false | randomly output one mapped position for ambiguous reads in a separated file |
| -u      | -unmapped | Boolean | false | output unmapped reads in a separated file |
| -A      | -ag-wild | Boolean | false | map using A/G bisulfite wildcards |
| -k      | -topk | Integer | 50 | maximum allowed mappings for a read in paried-end mapping|
| -L      | -fraglen | Integer | 1000 | max fragment length in paired-end mapping |

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

	walt -i hg19.dbindex -r read_1.fq -o reads_1_mapping.out
    
If mapping the reads from the *_2 reads file, the -A option should be set. This means that all Gs in the reads and genome are transfered to As. If -A option is not set, all Cs in the reads and genome are transfered to Ts.

    walt -i hg19.dbindex -r read_2.fq -A -o reads_2_mapping.out
    
Additionally, WALT supports comma-separated list of read files. WALT produces one mapping output file for each read file. For single-end mapping, the output file names will be appended "_s1", "_s2", and so on. Notice: except the first file path, all other file paths cannot be use '~'. For example, "-r ~/read_file1.fq,~/read_file2.fq" is not allowed. It should be "-r ~/read_file1.fq,/home/read_file2.fq", since linux system doesn't know it is a path except the first one.
	 
	 walt -i hg19.dbindex -r read_file1.fq,read_file2.fq,read_file3.fq -A -o reads_2_mapping.out

    
The default number of maximum allowed mismatches is 6. The maximum allowed mismatches can be set using -m option.

    walt -i hg19.dbindex -r read_1.fq -m 4 -o reads_1_mapping.out
    
The option -N sets the number of reads to mapping in each loop. If N is larger, the program takes large memory, especially for paired-end read mapping. If N is 1000000, both single-end and paired-end mapping take about 15 Gb memory. If N is 5000000, single-end mapping takes about 16 Gb memory, and paired-end mapping takes about 32 Gb memory. If N is set to be larger than 5000000, the program will set N to be 5000000 since when N is too large the program will take large memory but it will not be faster.
    
    walt -i hg19.dbindex -r read_1.fq -N 1000000 -o reads_1_mapping.out
    
To output the ambiguous mapped reads or unmapped reads, -u and -a options should be set. The ambigous mapped reads and unmapped reads will output to a separated file. For paired-end mapping, there will be a unmapped file and an ambiguous file for each of the mate 1 and mate 2 reads file.
    
    walt -i hg19.dbindex -r read_1.fq -u -a -o reads_1_mapping.out
    
For paired-end reads, -1 and -2 options are used for the mate reads files.
    
    walt -i hg19.dbindex -1 read_1.fq -2 read_2.fq -N 5000000 -o paired_reads_mapping.out
    
Similarly, WALT supports comma-separated list of read files for paired-end mapping. WALT produces one mapping output file for each read file pair. The output file names will be appended "_p1", "_p2", and so on. One other thing to note is mate 1 and mate 2 paired files should be in the same order.

	 walt -i hg19.dbindex -1 read_file1_1.fq,read_file2_1.fq,read_file3_1.fq \ 
	                      -2 read_file1_2.fq,read_file2_2.fq,read_file3_2.fq \
	                      -o paired_reads_mapping.out
    
    
    
### Contacts ###

***Haifeng Chen***   *haifengc@usc.edu*

***Andrew Smith***   *andrewds@usc.edu*

***Ting Chen***   *tingchen@usc.edu*
