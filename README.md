## bsmapper v1.0 ##

***bsmapper*** is a reads mapping program for bisulfite sequencing in DNA methylation study.


### Installation ###
(1) Download the source code from Github 

    git clone git@github.com:smithlabcode/bsmapper.git

(2) Build and Install
    
    cd bsmapper
    make all
    make install


### Indexing Genome ###
    
    makedb -c <genome folder or file> -o <index file>

### Bisulfite Mapping ###

single-end reads

    bsmapper -i <index file> -r <reads file> -o <output file> [options]

paired-end reads

    bsmapper -i <index file> -1 <reads file_1> -2 <reads file_2> -o <output file> [options]


### Mapping Options ###


| Option | Long Tag | Type | Default | Brief Description |
| :-------------: |:-------------:|:-----:|:-----:| :-----|
| -i      | -index | String | NULL |index file created by ***makedb*** command ( .dbindex) |
| -r      | -reads | String | NULL | single-end reads file (.fastq or .fq) |
| -1      | -reads1 | String | NULL | paired-end reads _1 file (.fastq or .fq) |
| -2      | -reads2 | String | NULL | paired-end reads _2 file (.fastq or .fq) |
| -o      | -output | String | NULL | output file name |
| -m      | -mismatch | Integer | 6 | maximum allowed mismatches |
| -N      | -number | Integer | 5000000 | number of reads to map at one loop |
| -A      | -ag-wild | Boolean | false | map using A/G bisulfite wildcards |
| -k      | -topk | Integer | 50 | maximum allowed mappings for a read in paried-end mapping|
| -L      | -fraglen | Integer | 1000 | max fragment length in paired-end mapping |

To see the list of options, use "-?" or "-help".

### Examples ###

(1) **Indexing Gneome**

For example, to make an index for UCSC hg19

	makedb -c hg19/ -o hg19.dbindex
   
or to make an index for chromosome 2

	makedb -c chr2.fa -o chr2.dbindex

The suffix of the index file should be '.dbindex'.
    
(2) **Bisulfite Mapping**

For example, to mapping reads to human genome hg19

	bsmapper -i hg19.dbindex -r read_1.fq -o reads_1_mapping.out
    
If mapping the reads from the *_2 reads file, the -A option should be set. This means that all Gs in the reads and genome are transfered to As. If -A option is not set, all Cs in the reads and genome are transfered to Ts.

    bsmapper -i hg19.dbindex -r read_2.fq -A -o reads_2_mapping.out
    
The default number of maximum allowed mismatches is 6. The maximum allowed mismatches can be set using -m option.

    bsmapper -i hg19.dbindex -r read_1.fq -m 4 -o reads_1_mapping.out
    
The option -N sets the number of reads to mapping in each loop. If N is larger, the program takes large memory, especially for paired-end read mapping. If N is 1000000, both single-end and paired-end mapping take about 15 Gb memory. If N is 5000000, single-end mapping takes about 16 Gb memory, and paired-end mapping takes about 32 Gb memory.
    
    bsmapper -i hg19.dbindex -r read_1.fq -N 1000000 -o reads_1_mapping.out
    
For paired-end reads, -1 and -2 options are used for the mate reads files.
    
    bsmapper -i hg19.dbindex -1 read_1.fq -2 read_2.fq -N 5000000 -o paired_reads_mapping.out
    
    
    
### Contacts ###

***Haifeng Chen***   *haifengc@usc.edu*


5/14/2015
