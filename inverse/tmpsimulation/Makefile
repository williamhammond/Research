OBJS = Main.o external.o matEg.o
CC = g++ -g
CFLAGS =-c -Wno-write-strings
OBJS = external.o matEg.o Main.o TMPSimulation.o Matrix.o
INC = -I/home/wth4280/lib -I/usr/local/include/matlab
LIB = -L/opt/MATLAB/R2015b/bin/glnxa64 -lgsl -lmx -leng -lgslcblas
Main : $(OBJS)
	    $(CC) $(OBJS) -o Main $(LIB) 
external.o : external.h
	    $(CC) $(CFLAGS) $(INC) external.cpp 
matEg.o : matEg.h inclu.h
	    $(CC) $(CFLAGS) $(INC) matEg.cpp 
TMPSimulation.o : TMPSimulation.h
	    $(CC) $(CFLAGS) $(INC) TMPSimulation.cpp 
Matrix.o: Matrix.h 
	    $(CC) $(CFLAGS) $(INC) Matrix.cpp 
Main.o : external.h inclu.h TMPSimulation.h matEg.h
	    $(CC) $(CFLAGS) $(INC) Main.cpp 

clean:
	    rm *.o Main
