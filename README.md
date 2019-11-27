MapCaller: An efficient and versatile approach for short-read mapping and variant identification using high-throughput sequenced data.
===================

Developers: Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu. Institute of Information Science, Academia Sinica, Taiwan.

# Introduction
MapCaller aligns every NGS short read against a reference genome and collects all the alignment information to deduce sequence variants. MapCaller adopts the mapping algorithm of KART to perform read alignment. It maintains a position frequency matrix to keep track of every nucleotide¡¦s frequency at each position in the reference genome, and it collects all insertion and deletion events which are found during the read mapping. Furthermore, MapCaller also learns all possible break points from discordant or partial read alignments. Finally, MapCaller identifies sequence variants based on all the above-mentioned information. The novelty of our algorithm derives from the integration of read mapping and the variation information gathering into a coherent system for genomic variant calling. Thus, the inputs to MapCaller is one or more NGS read file and an index file for the reference genome, and the output is a VCF file for the variant calling result.

For more information, please refer to our manuscript (https://www.biorxiv.org/content/10.1101/783605v2)

The benchmark data sets and the performance evaluation program could be found at http://bioapp.iis.sinica.edu.tw/~arith/MapCaller/

# Installation

## Conda
Install [Bioconda](https://bioconda.github.io/user/install.html) then type:
```
$ conda install -c conda-forge -c bioconda -c defaults mapcaller
```

## Homebrew
Install [HomeBrew](http://brew.sh/) (MacOS) or [LinuxBrew](http://linuxbrew.sh/) (Linux) then type:
```
$ brew install brewsci/bio/mapcaller
```

# Download

Please use the command 
  ```
  $ git clone https://github.com/hsinnan75/MapCaller.git
  ```
to download the package of MapCaller.

# Dependencies
To compile MapCaller, it requires `libboost-all-dev`, `libbz2-dev`, and `liblzma-dev` pre-installed in your system.

# Compiling
To compile MapCaller and the index tool, please just type `make` to compile MapCaller. If the compilation or the program fails, please contact me (arith@iis.sinica.edu.tw) or report it at the
[Issue Tracker](https://github.com/hsinnan75/MapCaller/issues).

# Test
You may run `run_test.sh` to test MapCaller with a toy example.

# Get updates
  ```
  $ bin/MapCaller update
  ```
or
  ```
  $ git fetch
  $ git merge origin/master master
  $ make
  ```

# Instructions

  ```
  $ bin/MapCaller [options]
  ```

# Usage

To index a reference genome, it requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

  ```
  $ bin/MapCaller index ref_file[ex.ecoli.fa] index_prefix[ex. Ecoli]
  ```
The above command is to index the genome file and store the index files begining with $index_prefix.

or

  ```
  $ bin/MapCaller -r ref_file[ex.ecoli.fa] 
  ```
The above command is to read a reference file and build a temporary index file

To perform variant calling, MapCaller requires the the index files of the reference genome and at least one read file (two read files for the separated paired-end reads). Users should use -i to specify the prefix of the index files (including the directory path).

 case 1: standard vcf output / sam output (optional) / bam output (optional)
  ```
 $ bin/MapCaller -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -vcf out.vcf [-sam out.sam][-bam out.bam]
  ```

 case 2: multiple input 
  ```
 $ bin/MapCaller -i ecoli -f ReadFileA_1.fq ReadFileB_1.fq ReadFileC_1.fq -f2 ReadFileA_2.fq ReadFileB_2.fq ReadFileC_2.fq -vcf out.vcf
  ```

 case 3: given a reference genome
  ```
 $ bin/MapCaller -r ecoli.fa -f ReadFile1.fa -f2 ReadFile2.fa -vcf out.vcf [-sam out.sam][-bam out.bam]
  ```

# File formats

- Reference genome files

    All reference genome files should be in FASTA format.

- Read files

    MapCaller is designed for NGS short reads. All reads files should be in FASTA/FASTQ format. Input files can be compressed with gzip format.
    Please note, if reads are stored in FASTA format, each read sequence is not allowed to be wrapped (split over multiple lines).
    Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. 
    If paired-end reads are separated into two files, use -f and -f2 to indicate the two filenames. The i-th reads in the two files are paired. If paired-end reads are in the same file, use -p. The first and second reads are paired, the third and fourth reads are paired, and so on. For the latter case, use -p to indicate the input file contains paired-end reads.

- Output file

    MapCaller outputs a SAM/BAM file [optional] and a VCF file. 

# Parameter setting

 ```
-t INT number of threads [16]

-i STR index prefix [optional, BWT based (BWA)]

-r STR reference filename [optional, fasta]

-f STR read filename [required, fasta or fastq or fq.gz]

-f2 STR read filename2 [optional, fasta or fastq or fq.gz], f and f2 are files with paired reads

-p the input read file consists of interleaved paired-end sequences [false]

-sam STR SAM output [optional, default: no mapping output]

-bam STR BAM output [optional, default: no mapping output]

-alg STR gapped alignment algorithm [optional, nw|ksw2, default: nw]

-vcf STR VCF output [output.vcf]

-dup INT Maximal PCR duplicates [optional, 1-100, default: 5]

-ploidy INT number of sets of chromosomes in a cell [optional, default:2]

-filter Apply variant filters (under test) [false]

-no_vcf No VCF output [false]

-gvcf GVCF mode [false]

-size Sequencing fragment size [default: 500, MapCaller can predict the fragment size automatically]

-ad INT Minimal ALT allele count [3]

-somatic detect somatic mutations [false]

-m output multiple alignments [false]

-v version number

-h help

  ```
# Changes
- 0.9.9.1: Adjusted read depth threshold.
- 0.9.9.2: Removed variants which appear in repetitive regions.
- 0.9.9.3: Adjusted read depth threshold for somatic mutation detection.
- 0.9.9.4: Added variant filters.
- 0.9.9.5: Fixed a bug on read count.
- 0.9.9.6: Fixed a bug on read count.
- 0.9.9.7: Fixed typos and warnings.
- 0.9.9.8: Used ksw2 algorithm for gapped alignments and fixed a bug on VCF format output.
- 0.9.9.9: Added a new argument(-ploidy) for setting different ploidy.
- 0.9.9.a: Fixed a bug on read mapping && Added an option (-dup) to filter out PCR-duplicates.
- 0.9.9.b: Added GVCF format (-gvcf).
- 0.9.9.c: Added a new argument(-monomorphic) to report all loci which do not have any potential alternates.
- 0.9.9.d: Fixed bugs on checking reference nucloetides and reporting deletion events.
- 0.9.9.e: Fixed bugs on reporting monomorphic regions and reported unmapped regions (potential deletions).
- 0.9.9.15: Fixed errors on VCF headers.
- 0.9.9.16: Added a new argument (-maxmm) for setting maximal mismatches in read alignment.
- 0.9.9.17: Updated VCF flags and format.
- 0.9.9.18: Fixed bugs && Added a new argument (-maxclip) for setting maximal clip size in read alignment.
- 0.9.9.19: Fixed a bug on the determination of deletion locus.
- 0.9.9.21: Reported CNV (N>=2) regions with the flag "DUP" (experimental).
- 0.9.9.22: Change the argument "-maxmm" to the maximal mismatch rate allowed in a read alignment.
- 0.9.9.23: Change the VCF output to fix the problem when applying vcftools with --plink
- 0.9.9.24: Summarize the process to a log file
- 0.9.9.25: Fixed typos in vcf headers.
- 0.9.9.26: Fixed a bug when running multi-threads on MacOS.
- 0.9.9.27: Integrated bwt_index into MapCaller
- 0.9.9.28: Build reference index on the fly

# Acknowledgements
We would like to thank [A/Prof. Torsten Seemann](https://github.com/tseemann)
and [Dr. Devon Ryan](https://github.com/dpryan79) for their help with the software.

# Note
The source code of the homemade genome sequence variation simulator is available at src/sv_simulator/.
