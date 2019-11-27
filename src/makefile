.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64 -msse4.1
LIB		= -lz -lm -lbz2 -llzma -lpthread -lstdc++
HTSLIB		= htslib
BWTLIB		= BWT_Index
SOURCE		= main.cpp GetData.cpp VariantCalling.cpp ReadMapping.cpp AlignmentRescue.cpp ReadAlignment.cpp AlignmentProfile.cpp SamReport.cpp tools.cpp bwt_index.cpp bwt_search.cpp nw_alignment.cpp ksw2_alignment.cpp KmerAnalysis.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o MapCaller $(BWTLIB)/libbwa.a $(HTSLIB)/libhts.a $(LIB) 

clean:
		rm -f *.o *~
