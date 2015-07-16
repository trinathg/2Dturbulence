DebugFlag=-pg -fopenmp -O3
ObjectFiles=driver.o initialize.o mesh.o tdma.o residual.o bcs.o solver.o visualize.o restart.o multigrid_v2p.o mulvec.o
Compile=g++
MKLPATH=
MKLINCLUDE=
MKLFLAG= 
Libraries=-larmadillo 
 
2Dturbsim: $(ObjectFiles)
	$(Compile) -o 2Dturbsim $(ObjectFiles) $(DebugFlag) $(MKLFLAG) $(Libraries)

driver.o: driver.cpp headers.h
	$(Compile) -c driver.cpp  $(DebugFlag) $(MKLFLAG)  

initialize.o: initialize.cpp headers.h init_2.h declarations.h
	$(Compile) -c initialize.cpp  $(DebugFlag) $(MKLFLAG)
 
mesh.o: mesh.cpp headers.h init_2.h declarations.h
	$(Compile) -c mesh.cpp  $(DebugFlag) $(MKLFLAG)

tdma.o: tdma.cpp headers.h init_2.h declarations.h
	$(Compile) -c tdma.cpp  $(DebugFlag) $(MKLFLAG)
 
residual.o: residual.cpp headers.h init_2.h declarations.h
	$(Compile) -c residual.cpp  $(DebugFlag) $(MKLFLAG)

bcs.o: bcs.cpp headers.h init_2.h declarations.h
	$(Compile) -c bcs.cpp  $(DebugFlag) $(MKLFLAG)

solver.o: solver.cpp headers.h init_2.h declarations.h
	$(Compile) -c solver.cpp  $(DebugFlag) $(MKLFLAG)

visualize.o: visualize.cpp headers.h init_2.h declarations.h
	$(Compile) -c visualize.cpp  $(DebugFlag) $(MKLFLAG)

restart.o: restart.cpp headers.h init_2.h declarations.h
	$(Compile) -c restart.cpp  $(DebugFlag) $(MKLFLAG)
	
multigrid_v2p.o: multigrid_v2p.cpp headers.h init_2.h declarations.h
	$(Compile) -c multigrid_v2p.cpp $(DebugFlag) $(MKLFLAG)	
	
mulvec.o: mulvec.cpp headers.h init_2.h declarations.h
	$(Compile) -c mulvec.cpp $(DebugFlag) $(MKLFLAG)	
		
clean:
	rm -f *.o data_* x-vel_* y-vel_* testdate
	rm -f *~
	touch *
