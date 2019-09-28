.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lstdc++
HTSLIB		= htslib
SOURCE		= main.cpp GetData.cpp VariantCalling.cpp ReadMapping.cpp AlignmentRescue.cpp ReadAlignment.cpp AlignmentProfile.cpp SamReport.cpp tools.cpp bwt_index.cpp bwt_search.cpp nw_alignment.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o MapCaller $(LIB) $(HTSLIB)/libhts.a -lbz2 -llzma

clean:
		rm -f *.o *~
