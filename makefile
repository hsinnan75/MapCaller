.KEEP_STAT:

all:		MapCaller bwt_index

MapCaller: htslib
		make -C src
		mv src/$@ bin/

htslib:
		make -C src/htslib

bwt_index:
		make -C src/BWT_Index
		mv src/BWT_Index/$@ bin/

clean:
		rm -f MapCaller bwt_index
		make clean -C src
		make clean -C src/htslib
		make clean -C src/BWT_Index

