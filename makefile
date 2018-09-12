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

all:		main index

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o MapCaller $(LIB)

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .

eva:		VarEva.cpp
		$(Compiler) $(FLAGS) VarEva.cpp -o var_eva

sim:		SVsim.cpp
		$(Compiler) $(FLAGS) SVsim.cpp -o SVsim                

clean:
		rm -f *.o *~
