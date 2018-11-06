.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lstdc++
SOURCE		= main.cpp GetData.cpp VariantCalling.cpp ReadMapping.cpp AlignmentRescue.cpp ReadAlignment.cpp AlignmentProfile.cpp SamReport.cpp tools.cpp bwt_index.cpp bwt_search.cpp nw_alignment.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o MapCaller $(LIB)

clean:
		rm -f *.o *~
