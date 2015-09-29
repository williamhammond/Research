//external.h
//declaration of external functions

#ifndef _EXTERNAL__H_
#define _EXTERNAL__H_

class Matrix;
//************************************************************************************************
//declaration of external functions
//************************************************************************************************
//---------------------------------------------------------------------------------------------
//test functions
//---------------------------------------------------------------------------------------------
/*
void anaTMP( char* filNam, double R, double* d,double** Um,int& inNum);

void anaBSP( char* filNam, double** Ut,int& outNum);

void anaTest ( char** filNams, int R, double* d);


//-------------------------------------------------------------------------------------------
//functions generating input files
//-------------------------------------------------------------------------------------------
void genStaFil( int staDim);

void genNoiFil(int staDim, int meaDim);
 */
// ‰»Î ‰≥ˆ∫Ø ˝
//ostream& operator<<(ostream&, FEM&);
//istream& operator>>(istream&, FEM&);
//------------------------------------------------------------------------------------
//Template function 
//--------------------------------------------------------------------------------------
//function name: swap

typedef struct{
	double x;
	double y;
	double z;
}cor;

void addPath(char** fullname,char* path, char** name, int num);
void writeFile(char* file, cor* data, int num, char* mode);
void wait();
void reSet(double* data, int len);

//void wriFile(char* file, mwArray& mat, char* mode);
//void test(char* filName,mwArray& mat, char* mode);
int inVol(cor a, cor b, cor c, cor d);
bool isoutTetra(cor gpCor,double* cor, int* tet, int num_tet);
double dotProduct(char axis, double x, double y, double z);
double dotProduct(double x,double y, double z, double x1, double y1, double z1);
void Poly3D(Matrix& p, double x, double y, double z);
void Poly3D_dx(Matrix& p, double x, double y, double z);
void Poly3D_dy(Matrix& p, double x, double y, double z);
void Poly3D_dz(Matrix& p, double x, double y, double z);
void testFECG();
int computeFECG();
void FECG(char** filNameBEM, char** filNamMFre, char** filOutput,char** TransOutput);
void FECG(char** filNameBEM, char** filNamMFree, char** filOutMFree, char** filOutBEM,char** TransOutput);


int checkNan(Matrix& p);
//-----------------------------------------------------------------
//applications
//-----------------------------------------------------------------
//void computeTorso(double*, double*,double*, int*, int, int, int, double*);
//------------------------------------------------------------------------------------------------------------------
//external function
//------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//d2c
//change coordinate in form of double array into cor
//-----------------------------------------------------------------------------------
inline cor d2c(double* d)
{
	cor c;
	c.x = d[0];
	c.y = d[1];
	c.z = d[2];
	return c;
}

inline void c2d(double* d, cor c)
{
	d[0] = c.x;
	d[1] = c.y;
	d[2] = c.z;
}


#endif