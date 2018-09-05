MapCaller: An efficient and versatile approach for short-read alignment and variant detection in high-throughput sequenced genomes
===================

Developers: Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu. Institute of Information Science, Academia Sinica, Taiwan.

# Introduction
MapCaller aligns every NGS short read against a reference genome and collects all the alignment information to deduce sequence variants. MapCaller adopts the mapping algorithm of KART to perform read alignment. It maintains a position frequency matrix to keep track of every nucleotide¡¦s frequency at each position in the reference genome, and it collects all insertion and deletion events which are found during the read mapping. Furthermore, MapCaller also learns all possible break points from discordant or partial read alignments. Finally, MapCaller identifies sequence variants based on all the above-mentioned information. The novelty of our algorithm derives from the integration of read mapping and the variation information gathering into a coherent system for genomic variant calling. Thus, the inputs to MapCaller is one or more NGS read file and an index file for the reference genome, and the output is a VCF file for the variant calling result.

The benchmark data sets and the performance evaluation program could be found at http://bioapp.iis.sinica.edu.tw/~arith/MapCaller/

# Download

Please use the command 
  ```
  $ git clone https://github.com/hsinnan75/MapCaller.git
  ```
to download the package of MapCaller.

# Compiling
To compile MapCaller and the index tool, please just type 'make' to compile MapCaller and bwt_index. If the compilation or the program fails, please contact me (arith@iis.sinica.edu.tw), Thanks.

# Changes
version 0.9.0: First version released.

# Get updates
  ```
  $ ./MapCaller update
  ```
or
  ```
  $ git fetch
  $ git merge origin/master master
  $ make
  ```

# Instructions

  ```
  $ ./MapCaller [options]
  ```

# Usage

To index a reference genome, it requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

  ```
  $ ./bwt_index ref_file[ex.ecoli.fa] index_prefix[ex. Ecoli]
  ```
The above command is to index the genome file Ecoli.fa and store the index files begining with ecoli.

Please note that if you find bwt_index does not work in your computer system, you may use bwa (http://bio-bwa.sourceforge.net/) to build the index files.
  ```
  $ ./bwa index -p index_prefix xxxx.fa
  ```
To perform variant calling, MapCaller requires the the index files of the reference genome and at least one read file (two read files for the separated paired-end reads). Users should use -i to specify the prefix of the index files (including the directory path).

 case 1: standard vcf output / sam output (optional)
  ```
 $ ./MapCaller -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -vcf out.vcf [-sam out.sam]
  ```

 case 2: multiple input 
  ```
 $ ./MapCaller -i ecoli -f ReadFileA_1.fq ReadFileB_1.fq ReadFileC_1.fq -f2 ReadFileA_2.fq ReadFileB_2.fq ReadFileC_2.fq -vcf out.vcf
  ```

# File formats

- Reference genome files

    All reference genome files should be in FASTA format.

- Read files

    All reads files should be in FASTA or FASTQ format. FASTQ files can be compressed with gzip format. We do not support FASTA files with gzip compression.
    Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. The alignment result will not be different in either format.

    If paired-end reads are separated into two files, use -f and -f2 to indicate the two filenames. The i-th reads in the two files are paired. If paired-end reads are in the same file, use -p. The first and second reads are paired, the third and fourth reads are paired, and so on. For the latter case, use -p to indicate the input file contains paired-end reads.

- Output file

    MapCaller outputs a SAM file [optional] and a VCF file. 

# Parameter setting

 ```
-t INT number of threads [16]

-i STR index prefix [BWT based (BWA), required]

-f STR read filename [required, fasta or fastq or fq.gz]

-f2 STR read filename2 [optional, fasta or fastq or fq.gz], f and f2 are files with paired reads

-p the input read file consists of interleaved paired-end sequences [default: false]

-sam STR SAM output [optional, default: no SAM output]

-vcf STR VCF output [default: output.vcf]

-no_vcf No VCF output [false]

-size Sequencing fragment size [default: 500, MapCaller can predict the fragment size automatically]

-m output multiple alignments [default: false]

  ```
