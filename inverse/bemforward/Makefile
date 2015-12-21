OBJS = Main.o external.o matEg.o
CC = g++ -g
CFLAGS =-c -Wno-write-strings
OBJS =  MFree.o Main.o matEg.o Matrix.o BEM.o external.o BEMCell.o 
INC = -I/usr/include/gsl/ -I/usr/local/include/matlab
LIB = -L/opt/MATLAB/R2015b/bin/glnxa64 -lgsl -lmx -leng -lgslcblas
Main : $(OBJS)
	    $(CC) $(OBJS) -o Main $(LIB) 
MFree.o: inclu.h BEMCell.h MFree.h BEM.h Matrix.h external.h
		$(CC) $(CFLAGS) $(INC) MFree.cpp
external.o : MFree.h external.h BEMCell.h inclu.h BEM.h matEg.h Matrix.h
	    $(CC) $(CFLAGS) $(INC) external.cpp 
matEg.o : matEg.h inclu.h Matrix.h
	    $(CC) $(CFLAGS) $(INC) matEg.cpp 
BEM.o : MFree.h inclu.h external.h BEMCell.h BEM.h
	    $(CC) $(CFLAGS) $(INC) BEM.cpp
Matrix.o: Matrix.h inclu.h external.h 
	    $(CC) $(CFLAGS) $(INC) Matrix.cpp 
BEMCell.o: inclu.h external.h BEMCell.h BEM.h
		$(CC) $(CFLAGS) $(INC) BEMCell.cpp
Main.o : external.h inclu.h BEM.h matEg.h
	    $(CC) $(CFLAGS) $(INC) Main.cpp 
clean:
	    rm *.o
		rm Main
