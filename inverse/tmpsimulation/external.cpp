#include "inclu.h"
#include "Matrix.h"
#include "external.h"

extern int ID;
extern int YDIM;
extern double DMAX;
extern int NUM_INF;
extern double BG_INTERVAL;
extern int IS_ISO;
extern int IS_HOMO;
extern char FHN_TYPE;
extern char PAR_TYPE;
extern int SCHEME;
extern char PATH_ROOT[MAX_FILENAME];

//***************************************
//check whether a matrix has nan
//****************************************
int checkNan(Matrix& p) {
	double* ptr = p.getDataPtr();
	int num = p.getCol() * p.getRow();
	for (int i = 0; i < num; i++) {
		if (isnan(ptr[i]))
			return 1;
	}
	return 0;
}

//**********************************************************
//Poly3D family for p, pdx...
//**********************************************************
void Poly3D(Matrix& p, double x, double y, double z) {
	p.setData(1.0, 0, 0);
	p.setData(x, 1, 0);
	p.setData(y, 2, 0);
	p.setData(z, 3, 0);
}

void Poly3D_dx(Matrix& p, double x, double y, double z) {
	p.setData(0.0, 0, 0);
	p.setData(1.0, 1, 0);
	p.setData(0.0, 2, 0);
	p.setData(0.0, 3, 0);
}

void Poly3D_dy(Matrix& p, double x, double y, double z) {
	p.setData(0.0, 0, 0);
	p.setData(0.0, 1, 0);
	p.setData(1.0, 2, 0);
	p.setData(0.0, 3, 0);
}

void Poly3D_dz(Matrix& p, double x, double y, double z) {
	p.setData(0.0, 0, 0);
	p.setData(0.0, 1, 0);
	p.setData(0.0, 2, 0);
	p.setData(1.0, 3, 0);
}

//-*******************************************************
//inVol
//called by "isout" 
//*******************************************************
int inVol(double* a, double* b, double* c, double* d) {
	double vol;
	double bxdx, bydy, bzdz, cxdx, cydy, czdz;

	bxdx = b[0] - d[0];
	bydy = b[1] - d[1];
	bzdz = b[2] - d[2];
	cxdx = c[0] - d[0];
	cydy = c[1] - d[1];
	czdz = c[2] - d[2];
	vol = (a[2] - d[2]) * (bxdx * cydy - bydy * cxdx)
			+ (a[1] - d[1]) * (bzdz * cxdx - bxdx * czdz)
			+ (a[0] - d[0]) * (bydy * czdz - bzdz * cydy);

	//cout << "vol: " << vol << endl;

	/* The volume should be an integer. */
	if (vol > 0)
		return 1;
	else if (vol < 0)
		return -1;
	else
		return 0;
}

//-*******************************************************
//inVol
//called by "isout" 
//*******************************************************
void reSet(double* data, int len) {
	for (int i = 0; i < len; i++)
		data[i] = 0;
}

void wait() {
	//cout<<"Press any key to continue"<<endl;
	//char c;
	//cin>>c;
	exit(0);
}

//*********************************************************
//reload function to compute cos(theta) between two vetors)
//*********************************************************
//ordinary
double dotProduct(double x, double y, double z, double x1, double y1,
		double z1) {
	return (x * x1 + y * y1 + z * z1)
			/ (sqrt(x * x + y * y + z * z) * sqrt(x1 * x1 + y1 * y1 + z1 * z1));
}

//with x,y or z axis, simplified
double dotProduct(char axis, double x, double y, double z) {
	switch (axis) {
	case 'x':
		return x / sqrt(x * x + y * y + z * z);
	case 'y':
		return y / sqrt(x * x + y * y + z * z);
	case 'z':
		return z / sqrt(x * x + y * y + z * z);
	default:
		cout << "error: nonexistent cartesian basis!" << endl;
		exit(0);
	}
}

//******************************
//round the double to the closest integer
//******************************
//int round(double d)
//{
//	return (int)floor(d);
//}

//*******************************
//read parameter
//********************************
void readPara() {
	printf("Enter the root path for input: \n");
//	scanf("%s\n",&input_path);
	gets(PATH_ROOT);

	cout << endl << "Input UDIM" << endl;
	cin >> ID;
	cout << "Input YDIM" << endl;
	cin >> YDIM;
	cout << endl << "Input DMAX" << endl;
	cin >> DMAX;
	cout << endl << "Input Number of nodes in influence domain" << endl;
	cin >> NUM_INF;
	cout << endl << "Input length for each background mesh" << endl;
	cin >> BG_INTERVAL;
	cout << endl << "Is Isotropic? (1, yes; 0, no)?" << endl;
	cin >> IS_ISO;
	cout << endl << "Is Homogeneous? (1, yes; 0, no)?" << endl;
	cin >> IS_HOMO;
	cout << endl << "TMP model scheme? 1, without stimuli; 2, with stimuli"
			<< endl;
	cin >> SCHEME;
	cout << endl << "FHN model type: P, M, O" << endl;
	cin >> FHN_TYPE;
	cout << endl << "Input parameter: a,others" << endl;
	cin >> PAR_TYPE;

	printf("*******check***********");
	printf("input root path = %s\n", PATH_ROOT);
	printf("UDIM = %d\n", ID);
	printf("YDIM = %d\n", YDIM);
	printf("DMAX = %f\n", DMAX);
	printf("NUM_INF = %d\n", NUM_INF);
	printf("BG_INTERVAL = %F\n", BG_INTERVAL);
	printf("IS_ISO = %d\n", IS_ISO);
	printf("IS_HOMO = %d\n", IS_HOMO);
	printf("SCHEME = %d\n", SCHEME);
	printf("MODEL TYPE = %c\n", FHN_TYPE);
	printf("PARAMETER = %c\n", PAR_TYPE);
}
