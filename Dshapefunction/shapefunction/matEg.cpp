/*
 * matEg.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: jxx2144
 */


#include "inclu.h"
#include "matEg.h"
#include "Matrix.h"



//********************************
//ctor/dtor
//********************************
matEg::matEg(int lenBuffer):nLenBuffer(lenBuffer)
{
	pBuffer = new char[nLenBuffer];
}

matEg::~matEg()
{
	delete []pBuffer;
}

//***************************************
//start the engine, begin the output buffer
//***************************************
int matEg::Open(const char* cmd)
{
	if (!(pEg = mat::engOpen(cmd)))
	{
		cout<<endl<<"Error: can not open the matlab engine!"<<endl;
		return EXIT_FAILURE;
	}

	cout<<endl<<"***Began matlab engine***"<<endl;
	engOutputBuffer(pEg,pBuffer,nLenBuffer);

	return 1;
}

//****************************************
//stop the engine, stop the output buffer
//****************************************
int matEg::Close()
{
	if(engClose(pEg))
	{
		cout<<endl<<"Error: can not stop the matlab engine or the engine has been closed!"<<endl;
		return EXIT_FAILURE;
	}

	engOutputBuffer(pEg,NULL,0);

	cout<<endl<<"***End matlab engine***"<<endl;

	return 1;
}

//****************************************
//for unique matrix in the whole program, such as sys.M
//from array needs to be computed in matlab, formulate mxArray, and put into workspace
//****************************************
//template <typename T>
int matEg::addData(double* data, const char* name,int row, int col)
{
/*	if (nCurrentMx >= nLenMx)
	{
		cout<<endl<<"Error: Matrix exceeds index!"<<endl;
		return EXIT_FAILURE;
	}
	else
	{
*/

	//creat matrix***********
	mat::mxArray* pMx = mat::mxCreateDoubleMatrix(row,col,mat::mxREAL);
	if(!pMx)
	{
		cout<<endl<<"Error: can not creat mxArray: "<<name<<endl;
		return EXIT_FAILURE;
	}
/*		//record its name***********
		pName[nCurrentMx] = new char[sizeof(name)];
		strcpy(pName[nCurrentMx],name);
*/		//copy the data***********
	memcpy((void*)mat::mxGetPr(pMx),(void*) data,row*col*sizeof(double));
	//put into matlab workspace*******
	if (engPutVariable(pEg,name,pMx))
	{
		cout<<endl<<"Error: can not put mxArray: "<<name<<" into matlab workspace"<<endl;
		return EXIT_FAILURE;
	}

	mat::mxDestroyArray(pMx);

	return 1;
}

//combined with Matrix**********************
int matEg::addData(Matrix& data, const char* name)
{
/*	if (nCurrentMx >= nLenMx)
	{
		cout<<endl<<"Error: Matrix exceeds index!"<<endl;
		return EXIT_FAILURE;
	}
	else
	{
*/

	int row = data.getRow();
	int col = data.getCol();
	double* ptr = data.getDataPtr();
	//creat matrix***********
	mat::mxArray* pMx = mat::mxCreateDoubleMatrix(row,col,mat::mxREAL);
	if(!pMx)
	{
		cout<<endl<<"Error: can not creat mxArray: "<<name<<endl;
		return EXIT_FAILURE;
	}
/*		//record its name***********
		pName[nCurrentMx] = new char[sizeof(name)];
		strcpy(pName[nCurrentMx],name);
*/		//copy the data***********
	memcpy((void*)mat::mxGetPr(pMx),(void*) ptr,row*col*sizeof(double));
	//put into matlab workspace*******
	if (engPutVariable(pEg,name,pMx))
	{
		cout<<endl<<"Error: can not put mxArray: "<<name<<" into matlab workspace"<<endl;
		return EXIT_FAILURE;
	}

	mat::mxDestroyArray(pMx);

	return 1;
}
//********************************************
//for the same matrice but related to a group of objects
//eg. shape for all gpts
//if the first time access, the same action as "addData"
//otherwise, change the value of the matrix
//*******************************************
/*
int matEg::changeData(double* data, const char* name,int row, int col)
{
	for (int i=0; i<nLenMx; i++)
	{
		if (strcmp(pName[i],name) == 0)
		{
			double* ptr = mxGetPr(pMx[i]);
			for(j=0; j<row*col; j++)
				ptr[j] = data[j];
            //put into matlab workspace
			if (engPutVariable(pEg,name,pMx[i]))
			{
				cout<<endl<<"Error: can not put mxArray #"<<i<<" "<<name<<" into matlab workspace"<<endl;
				return EXIT_FAILURE;
			}
		}
		else
			addData(data,name,row,col);
	return 1;
}
*/
//********************************************
//get the data from the result of matlab session, by the name
//********************************************
//template <typename T>
int matEg::getData(double* data, const char* name)
{
/*	for (int i=0; i<nLenMx; i++)
	{
		if (strcmp(pName[i],name) == 0)
		{
			//destroy early mxArray********important!!***
			if(pMx[i])
				mxDestroyArray(pMx[i]);
*/
	mat::mxArray* pMx = mat::engGetVariable(pEg,name);
	if (pMx  == NULL)
	{
		cout<<"Errors: can not read the variable"<<name<<"from matlab"<<endl;
		return EXIT_FAILURE;
	}
	double* ptr = mat::mxGetPr(pMx);
	int num = mat::mxGetM(pMx)*mat::mxGetN(pMx);
	for(int j=0; j<num; j++)
		data[j] = ptr[j];

	mat::mxDestroyArray(pMx);
	return 1;
}

//***combined with Matrix**********************
int matEg::getData(Matrix& data, const char* name)
{
/*	for (int i=0; i<nLenMx; i++)
	{
		if (strcmp(pName[i],name) == 0)
		{
			//destroy early mxArray********important!!***
			if(pMx[i])
				mxDestroyArray(pMx[i]);
*/
	mat::mxArray* pMx = mat::engGetVariable(pEg,name);
	if (pMx  == NULL)
	{
		cout<<"Errors: can not read the variable"<<name<<"from matlab"<<endl;
		return EXIT_FAILURE;
	}
	double* ptr = mat::mxGetPr(pMx);
	data.setData(ptr);
	mat::mxDestroyArray(pMx);
	return 1;
}

//*******************************************
//for a temptoray scale in Mat workspace, get double in C without cretaing matrix in engine
//******************************************
double matEg::getTemptScale(const char* name)
{
   double* ptr = mat::mxGetPr(mat::engGetVariable(pEg,name));
   return ptr[0];
}


//*************************************
//send command to matlab for computing
//*************************************
int matEg::evalString(const char* string)
{
	if (mat::engEvalString(pEg,string))
	{
		cout<<endl<<"Error: can not evaluate the function"<<endl;
		return EXIT_FAILURE;
	}

	return 1;
}





