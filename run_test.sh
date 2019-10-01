#!/bin/bash
#test indexing
echo
echo "Test1 -- Generate index files with a reference file"
echo "Command=bin/bwt_index test/ref.fa test/RefIdx"
echo
bin/bwt_index test/ref.fa test/RefIdx

#test alignment
echo
echo "Test2 -- Align paired-end reads with 4 threads"
echo "Command=bin/MapCaller -i test/RefIdx -t 4 -f test/r1.fq -f2 test/r2.fq -vcf test/out.vcf"
echo
bin/MapCaller -i test/RefIdx -t 4 -f test/r1.fq -f2 test/r2.fq -vcf test/out.vcf

echo
echo "[End of test]"
