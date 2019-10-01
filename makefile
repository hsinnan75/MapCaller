.KEEP_STAT:

all:		MapCaller bwt_index

MapCaller: htslib
		$(MAKE) -C src
		mkdir -p bin/ && cp -f src/$@ bin/

htslib:
		$(MAKE) -C src/htslib libhts.a

bwt_index:
		$(MAKE) -C src/BWT_Index
		mkdir -p bin/ && cp -f src/BWT_Index/$@ bin/

clean:
		rm -f MapCaller bwt_index
		$(MAKE) clean -C src
		$(MAKE) clean -C src/htslib
		$(MAKE) clean -C src/BWT_Index

