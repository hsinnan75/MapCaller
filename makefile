.KEEP_STAT:

all:		MapCaller

MapCaller: index htslib
		$(MAKE) -C src
		mkdir -p bin/ && cp -f src/$@ bin/
index:
		$(MAKE) -C src/BWT_Index libbwa.a
htslib:
		$(MAKE) -C src/htslib libhts.a

clean:
		rm -f bin/MapCaller
		$(MAKE) clean -C src
		$(MAKE) clean -C src/htslib
		$(MAKE) clean -C src/BWT_Index

