.KEEP_STAT:

all:		MapCaller bwt_index

MapCaller: htslib
		make -C src
		mkdir -p bin/ && mv src/$@ bin/

htslib:
		make -C src/htslib

bwt_index:
		make -C src/BWT_Index
		mkdir -p bin/ && mv src/BWT_Index/$@ bin/

clean:
		rm -f MapCaller bwt_index
		make clean -C src
		make clean -C src/htslib
		make clean -C src/BWT_Index

