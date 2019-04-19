.KEEP_STAT:

all:		htslib main index

htslib:
		make -C src/htslib

main:
		make -C src && mv src/MapCaller .

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .

