#ifndef EXTERNAL_H_
#define EXTERNAL_H_

class Matrix;


int checkNan(Matrix& p);

void Poly3D(Matrix& p, double x, double y, double z);
void Poly3D_dx(Matrix& p, double x, double y, double z);
void Poly3D_dy(Matrix& p, double x, double y, double z);
void Poly3D_dz(Matrix& p, double x, double y, double z);

int inVol(double* a, double* b, double* c, double* d);

void reSet(double* data, int len);
void wait();

double dotProduct(char axis, double x, double y, double z);
double dotProduct(double x,double y, double z, double x1, double y1, double z1);
//int round(double d);

void readPara();

#endif /*EXTERNAL_H_*/
