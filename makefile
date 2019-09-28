.KEEP_STAT:

all:		MapCaller bwt_index

MapCaller: htslib
		make -C src
		mv src/$@ .

htslib:
		make -C src/htslib

bwt_index:
		make -C BWT_Index
		mv BWT_Index/$@ .

clean:
		rm -f MapCaller bwt_index
		make clean -C src
		make clean -C src/htslib
		make clean -C BWT_Index

