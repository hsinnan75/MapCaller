CC=		gcc
FLAGS=		-I${PREFIX}/include -L${PREFIX}/lib -g -Wall -O3
AR=		ar
LOBJS=		utils.o bwt.o bntseq.o QSufSort.o bwt_gen.o bwtindex.o
PROG=		bwt_index
LIBS=		-lm -lz

.c.o:
		$(CC) -c $(FLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f *.o $(PROG) *~ *.a
