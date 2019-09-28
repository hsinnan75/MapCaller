.KEEP_STAT:

all:		htslib index main

htslib:
		make -C src/htslib

main:
		make -C src && mv src/MapCaller .

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .

clean:
		rm -f MapCaller bwt_index
		make clean -C src
		make clean -C src/htslib
		make clean -C BWT_Index

