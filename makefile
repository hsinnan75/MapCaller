.KEEP_STAT:

all:		main index

main:
		make -C src && mv src/MapCaller .

index:
		make -C BWT_Index && mv BWT_Index/bwt_index .

eva:		VarEva.cpp
		$(Compiler) $(FLAGS) VarEva.cpp -o var_eva

sim:		SVsim.cpp
		$(Compiler) $(FLAGS) SVsim.cpp -o SVsim                
