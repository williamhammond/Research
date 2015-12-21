#ifndef INCLU_H_
#define INCLU_H_

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
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
using std::cerr;
using std::cout;
using std::cin;
using std::endl;


//**************************************************************
//constant
//**************************************************************
//const int ID = 1721;//1230;//1373;//1239;//
const int DIM = 3;
//const double DMAX = 1.5;
//const double BG_INTERVAL =5;
const int MAX_NAT_GP = 10;
const int NUM_GP = 3;
const int MAX_GP = 70000;
//const int NUM_INF = 20;
const double MAX_R = 1000;
const int MAX_INF = 1000;
const int MAX_FILENAME = 500;
const int MAX_BUF = 5000;
//const int IS_ISO = 0;
const double TOL_BG = 1;     // RELAX THE BG MESH BY SOME MARGIN
const double PI = 3.141592657;

//note, in this program, D0 initialize as sqrt(D0)
const double DI = 1.0;
const double DI_L = 2.0;      // this is the square root of diffusion parameter
const double DI_S = DI_L/2.5819;   //square root of the anisotropic ratio of diffusion, i.e., the anistropic ration of CV
const double Ransio = 2.5819;
const double P2M = 4.0;    // ratio of conduction velocity of PJN to Myocardium.
const double DET_INVERSE = 1e-16;

const double U_INIT = 1.0;
const double U_THRESH_PEAK = 0.7;
const int INTERVAL_AT = 10;
const double THRESH = 0.3;
const double FOCI=0.5;
const int STI_DUR = 54;// 30; //   // stimulus duration, should be longer then RV delay
const int INTERVAL_RT_F = 4;
const int INTERVAL_RT_B =30;
const double THRESH_UP = 0.5;         //used to decide the t_s, t_f for detecting AT

const double PAR_AO = 0.1;
const double PAR_BO = 0.001;
const double PAR_DO = 3;
const double PAR_C1 = 0.26;
const double PAR_C2 = 0.1;
const double PAR_AM = 0.13;
const double PAR_BM = 0.013;
const double PAR_DM = 1.0;
const double PAR_K = 8;
const double PAR_E = 0.01;
const double PAR_AP = 0.15;



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
	strcpy(testFile,"Test/");
	strcat(testFile,filName);

	writeFile(testFile,data,size,num,mode);
}




#endif /*INCLU_H_*/
