
#ifndef MATEG_H_
#define MATEG_H_

//do not take any memories except buffer, only used as a tool to put double* or Matrix into MATLAB and getbak
//#include <stdlib.h>
//#include <stdio.h>
//#include "inclu.h"
//************************************************
#include "engine.h"

using std::ostream;
using std::cout;
using std::cerr;
using std::endl;
class Matrix;

//************************************************
class matEg
{
public:
	matEg(int lenBuffer);
	~matEg();

	int Open(const char*);
	int Close();
//	template<typename T> 
	int addData(double* data, const char* name,int row, int col);
	int addData(Matrix& m, const char* name);
//	int changeData(double* data,const char* name);
//	template<typename T> 
	int getData(double* data, const char* name);
	int getData(Matrix& m, const char* name);
	double getTemptScale(const char* name);
	int evalString(const char *string);
	  
private:
	Engine* pEg;
//	mxArray** pMx;                      //matlab matrix used to transfer data between C & matlab
//   char** pName;                   // corresponding name of these matrix
//	int nLenMx;
//	int nCurrentMx;
    char* pBuffer;     // one in each engine, start at the beginning, close at the end
	int nLenBuffer;
};


#endif /*MATEG_H_*/
