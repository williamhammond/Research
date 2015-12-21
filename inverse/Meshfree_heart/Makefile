OBJS = Main.o external.o matEg.o
CC = g++ -g
CFLAGS =-c -Wno-write-strings
OBJS = external.o Heart.o matEg.o Matrix.o main.o
INC = -I/usr/include/gsl/ -I/usr/local/include/matlab
LIB = -L/opt/MATLAB/R2015b/bin/glnxa64 -lgsl -lmx -leng -lgslcblas

Main : $(OBJS)
	    $(CC) $(OBJS) -o Main $(LIB) 
Heart.o : inclu.h external.h Heart.h matEg.h Matrix.h
		$(CC) $(CFLAGS) $(INC) Heart.cpp
external.o : external.h Matrix.h Heart.h inclu.h
	    $(CC) $(CFLAGS) $(INC) external.cpp 
matEg.o : matEg.h inclu.h Matrix.h
	    $(CC) $(CFLAGS) $(INC) matEg.cpp 
Matrix.o: inclu.h external.h Matrix.h
		$(CC) $(CFLAGS) $(INC) Matrix.cpp
main.o : external.h inclu.h matEg.h
	    $(CC) $(CFLAGS) $(INC) main.cpp 

clean:
	    rm *.o main
