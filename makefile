.KEEP_STAT:

all:		MapCaller bwt_index

MapCaller: htslib
		$(MAKE) -C src
		mkdir -p bin/ && cp -f src/$@ bin/

htslib:
		$(MAKE) -C src/htslib libhts.a

clean:
		rm -f bin/MapCaller
		$(MAKE) clean -C src
		$(MAKE) clean -C src/htslib

