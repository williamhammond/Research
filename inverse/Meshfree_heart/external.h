/*
 * external.h
 *
 *  Created on: Jul 3, 2012
 *      Author: jxx2144
 */

#ifndef EXTERNAL_H_
#define EXTERNAL_H_

class Heart3D;
class Heart;
class Matrix;
class U;
//*****************************************************************************
// functions
//*****************************************************************************
//***************************************
//read files
//**************************************

template<typename T>
void readFile(char* file, T** data, size_t size, int col_dim,int& row_dim)
{
	//char* file = strcpy(dirNam,filNam);
	FILE* hfil = fopen(file,"rb");
	if(hfil==NULL)
	{
		::cerr<<" File open failure!"<<file<<endl;
		exit(0);
	}
	else
	{
		fseek(hfil,0,SEEK_END);
        int filSiz = ftell(hfil);
		rewind(hfil);
        //get the size(rowDim*colDim)of input matrix
		int datByt = filSiz/size;
		row_dim = datByt/col_dim;
		(*data) = new T[datByt];
		fread((*data),size,datByt,hfil); // read out colDim data for every row, each data size is specified by"size"
	}
	fclose(hfil);
}

//*******************************************
//write files
//*******************************************
template<typename T>
void writeFile( char* file, T data, size_t size, int num, char* mode)
{
	FILE* hfil = fopen(file,mode);
	fwrite(data,size,num,hfil);
	fclose(hfil);
}

//*******************************************
//write testing files
//*******************************************
/*template<typename T>
void test(char* file,T data, size_t size, int num,char* mode)
{
	char* testFile = new char[MAX_FILENAME];
	strcpy(testFile,"../3DHeartModeling/test/");
	strcat(testFile,file);

	writeFile(testFile,data,size,num,mode);
}*/

//****************************************
//set a double array to be all 0s
//****************************************
template<typename T>
void reSet(T* data, int num)
{
	int i = 0;
	for(i=0; i<num; i++)
		data[i] = 0;
}

//*************************************
// tools
//**************************************
/*char* addPath(char* path, char* name);

void  poly3D(Matrix& p,double x, double y, double z);
void  poly3D_dx(Matrix& p,double x, double y, double z);
void poly3D_dy(Matrix& p,double x, double y, double z);
void poly3D_dz(Matrix& p,double x, double y, double z);

*/
/*
double Double(const mwNumericSubArray A);

double Double(const mwArray A);

mwArray operator~ (const mwArray A);

mwArray operator += (mwArray& A, const mwArray B);

mwArray operator -= (mwArray& A, const mwArray B);

mwArray operator *= (mwArray& A, const mwArray B);
*/
//void Newmark(NewM& nm,double delta_t);
void get_uI(U *uI_t, U &uI_k);

//**********************************
//inVol: test whether a point is in a tetra
//       used in 2DHeart::InTetra
//**********************************
int inVol(tPointd a, tPointd b, tPointd c, tPointd d);
//**********************************
//whether point pt in the tetras defined by cor, tet with num elements
//**********************************
bool outTetra(tPointd pt, double* cor, int* tet, int num);

//**********************************
//dotproduct to compute cos between two vectors, which share c as common origin
//***********************************
double dotProduct(double ax, double ay, double az,
				  double bx, double by, double bz,
				  double cx, double cy, double cz);

double dotProduct(tPointd r1, double* r2);
//**********************************
//crossproduct to compute sin between two vectors, which share c as common origin
double crossProduct(double ax, double ay, double az,
				  double bx, double by, double bz,
				  double cx, double cy, double cz,
				  double ox,double oy, double oz);
void crossProduct(tPointd r, tPointd r1,tPointd r2);

//**************************************
//compute distance in 3D space between pts r0 and line r1-r2
//***************************************
double distancePtsLine3D(tPointd r0, tPointd r1, tPointd r2);

//*******************************
//norm of vector r
//*******************************
double normVec(tPointd r);

//******************************
//arc tan in 2PI range
//******************************
double atan2PI(double x, double y);

//**************************************
//inverse of 3*3 matrix
//**************************************
void invMatrix3(double* inv, double* data);
void invMatrix3(double* data);

//***************************************
//generate mfree particles given heart surface geometry
//***************************************
void genMFreeNodes(char*,char*);

//*****************************************
//initial processing
//*****************************************
void preProcessing(char* file_in, char* file_out,char*file_image);

//*******************************************
//fiber mapping given reference heart geometry and fiber structure
//                    & target heart geometry
//********************************************
void fibMapping(char* file_in_source, char* file_in_target, char* file_out,char* file_image);

//*******************************************
//Another way of fiber mapping: surface deform + interp
//********************************************
void fibMappingInterp(char* file_in_source, char* file_in_target, char* file_out);

//********************************************
//from ref, get the function between tdr & tdtheta
//********************************************
void fibTrainning(char* file_in, char* file_out);

//********************************************
//find the source points in target heart, given those in ref heart
//*******************************************
void finalModeling(char* file_ref,char* file_new,int fib_ready);
//********************************************
//given rotated heart geo & fib, generate fib in original cor system
//*******************************************
void getFinalFiber(char*, char*, char*);

//*****************************************
//fibReMapping: based on reigistered heart, fiber mapping by epi/l/rendo/mfree
//*****************************************
void fibReMapping(char* file_in, char* file_ref, char* file_out);
void getClosestPts(int& seq_out, double& r_out, int seq_in, double* in_cor, int* ref_idx, int ref_num, int type);
void determineP1P2(double& dr, double& tr, int& id_P1, int& id_P2, double r_epi, int id_epi, double r_lendo, int id_lendo, double r_rendo, int id_rendo, int id, double* cor);

void rigidReg(double* cor, int mn_num_nodes);

//****************************************
//surfaceSmoothing: use 2-pass Gaussian filter
//***************************************
void surfaceSmoothing(char* path,double lamda1, double lamda2, int n);
//****************************************
//surfaceSmoothing:change original lv_tet to e useable in smoothed data
//***************************************
void getSmoothLV(char* path, bool addlendo);
//****************************************
//rotate the geometry
//***************************************
void geoRotate(char* path_in, char* path_ref);

void dataSetRotate(double* cor,int num, tPointd rotAng);
//***************************************
//tools for fibReMappng, group heart nodes
//***************************************
void groupHeartNodes(double** l, int** sl, int nl, double* all, int n_all, int* idx, int type);  //all nodes

//***************************************
//tools for fibReMapping, given each subset of nodes, find the cloeset points,
//and idx them back to original all nodes
//****************************************
void findClosestPart(int* seq_all,
		             double** cor_target, double** cor_ref,
		             int** seq_target,int** seq_ref,
		             int num_target, int num_ref,
		             double* all_in, double* all_ref,
		             int n_in, int n_ref,
		             int* idx_in, int* idx_ref,
		             int type);

//*****************************************
//find closest points
//*****************************************
void getClosestPts(double** closest, double* cor_ref,double* cor_target, int num_ref, int num_target);

//*********************************************
//getNumOfSubsets: get the number of lendo,rendo, epi & mfree
//*********************************************
void getNumOfSubsets(int& l,int& r, int& e, int& m, int* idx, int n);

//*****************************************
//show the result of final modeling
//****************************************
void showResult(Heart& heart);

//******************************
//given a ref tetra with surface labelled & a meshfree data with seg labelled, get the depth division
//*****************************
void transmuralDivision(char* file_ref, char* file_new);

//*******************************************
//visualization part (vtk)
//*******************************************
/*void pointSetConstruction(vtkPolyData* poly,double* cor, int dim_c);

void getData(double* cor,vtkDataSet* dataSet, int dim);

void dataSetUpdate(vtkPolyData* poly, double* cor, int num);*/
//****************************************
//other tools
//****************************************
//void reSet(double* data, int num);

void wait();






#endif /* EXTERNAL_H_ */
