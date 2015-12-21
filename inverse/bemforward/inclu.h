#ifndef _INCLU__H_
#define _INCLU__H_

//**************************************************************
//inclusions
//**************************************************************
#include<iostream>
#include<iomanip>
#include<string>
#include <fstream>
#include<stdlib.h>
#include<new>
#include<time.h>
#include <math.h>
#include<cstring>
using std::cerr;
using std::cout;
using std::cin;
using std::endl;

//**************************************************************
//constant
//**************************************************************
//for BEM&MFree********************
//control************
const int NUM_NODE = 1898;//
const int GAUPOI = 7;//
const int SINGAUPOI = 3;
const int GPNUM = 3;
const int INFNUM = 15;//30; //
const double BG_INTERVAL =6;//6;//0.15;//
const double DMAX = 1.2;
const int FLAG = 0;                    //whether class MFree generate B again ( 0 == yes)
const bool ISO = false;     // whether iso (=1) / aniso (=0)
const bool FORBIDDEN= false;
const bool FIBER_INTER = true;  //gp fib interpolated from meshfree nodes
const bool FIBER_CLOSE = false; //gp fib interpolated from closest nodes
const int IS_SIN = 1;
const char OUT_BOUNDARY = 'T';
const int IS_CONVEX = 1;
const double DET_INVERSE = 0.00001;
//memeory*********************
const int MAX_GP = 80000;
const int MAX_INF = 300;
const int MAX_NAT_GP = 10;
const double MAX_R = 100000;
const int MAX_FILENAME = 100;
//const double MAX_COR = 0.5;
const int MAX_BUF = 5000;
const double PI =  3.1415926;
const int DIM = 3;
const double NORM_SCALE = 10;
const double SCALE_T = 2;
//const double Dl = 4.0;
//const double RT = 4.0;
//***material property****************
const double RIe = 0.3;
const double RKe = 0.8;  //ratio for approximating isotropic bulk
const double DT = 0.2*0.001;  // con 
const double DK_L = 0.48*0.001;
const double DK_S = 0.12*0.001;
const double DI_L = 0.24*0.001;
const double DI_S = 0.024*0.001;
const double DI = RIe*DI_L + (1-RIe)*DI_S;
const double DK = RKe*DK_L + (1-RKe)*DK_S;           // real conductivity value
/*
const double Deff = DK;                  // used in F2
const double DTe = DT;
const double DKe = DK;
*/
const double Deff = (DK+DT)/2.0;//(5*DK + 7*DT)/12.0;// (10*DK - 7*DT)/3.0;          // F1
const double DTe = Deff;
const double DKe = Deff;
//**file********************************
extern char* FILE_REF_COR; 
extern char*  FILE_REF_FIB; 

//**************************************************************
//tmplate functions
//**************************************************************
//**numeric tools*****************
template<class T>
void swap(T& t1, T& t2)
{
	T temp;
	temp = t1;
	t1 = t2;
	t2 = temp;
}

//-*****************
//function name: absolute
//*****************
template<class T>
T absolute(T t)
{
	if (t<0)
		t = -t;
	return t;
}

//*****************
template<class T>
T power(T t)
{
	return t*t;
}

//file tools********************************
template<typename T>
void readFile(const char* filNam, T** data, size_t size,int colDim,int& rowDim)
{
	if(*data)
		delete [](*data);
	
	FILE* hfil = fopen(filNam,"rb");
	if(hfil==NULL)
	{
		::cerr<<" File open failure!"<<filNam<<endl;
//		wait();
		exit(0);
	}
	else
	{
		fseek(hfil,0,SEEK_END);
        int filSiz = ftell(hfil);
		rewind(hfil);
        //get the size(rowDim*colDim)of input matrix
		int datByt = filSiz/size;
		rowDim = datByt/colDim;
		(*data) = new T[datByt];
		fread((*data),size,datByt,hfil); // read out colDim data for every row, each data size is specified by"size"
	}
	fclose(hfil);
}
    
//*****************
template<typename T>
void writeFile(const char* filNam, T data, size_t size, int num, char* mode)
{	
	FILE* hfil = fopen(filNam,mode);
	fwrite(data,size,num,hfil);
	fclose(hfil);
}

//*****************
template<typename T>
void test(const char* filName,T data, size_t size, int num,char* mode)
{
	char* testFile = new char[MAX_FILENAME];
	strcpy(testFile,"../BEMForward/Test/");
	strcat(testFile,filName);

	writeFile(testFile,data,size,num,mode);
}


#endif //_INCLU__H_
