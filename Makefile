# A Makefile for bayescor, bbnet and gbnet 

objcomm = bayesub.o globals.o sa.o CmdLine.o fisher2.o
objfunc = func.o CmdLine.o 
objscor = bayescor.o $(objcomm)
objbb = bbnet.o $(objcomm)
objgb = gbnet.o $(objcomm)

func bayescor bbnet gbnet: $(objfunc) $(objscor) $(objbb) $(objgb)
	g++ -m32 -o func $(objfunc)
	g++ -m32 -o bayescor $(objscor)
	g++ -m32 -o bbnet $(objbb)
	g++ -m32 -o gbnet $(objgb)

func.o: prepsub.h CmdLine.h
	g++ -O3 -m32 -c func.cpp
bayescor.o bbnet.o: bayesub.h globals.h CmdLine.h
	g++ -O3 -m32 -c bayescor.cpp bbnet.cpp
gbnet.o: bayesub.h globals.h sa.h CmdLine.h
	g++ -O3 -m32 -c gbnet.cpp
bayesub.o: bayesub.h globals.h sa.h fisher2.h
	g++ -O3 -m32 -c bayesub.cpp
sa.o: sa.h badefs.h
	g++ -O3 -m32 -c sa.cpp
globals.o: globals.h 
	g++ -O3 -m32 -c globals.cpp
CmdLine.o: CmdLine.h
	g++ -O3 -m32 -c CmdLine.cpp
fisher2.o: fisher2.h Boolean.h Constants.h Memory.h
	g++ -O3 -m32 -c fisher2.cpp

clean:
	rm -f $(objfunc) $(objscor) $(objbb) $(objgb)
	rm -f func bayescor bbnet gbnet

