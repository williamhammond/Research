#ifndef MATRIX_H_
#define MATRIX_H_

/******************************************************************
 //contain most basic matrix operations
 ******************************************************************/

class Matrix {
	//-public member data
public:
	//-constructor, copy constructor, and deconstructor
	Matrix();
	Matrix(double* data, int row, int col);
	Matrix(int row, int col);
	Matrix(const Matrix& m);
	Matrix(char* file, int dim, char roc);
	~Matrix();

	//-Matrix operations.
	//    double det();
	void operator=(const Matrix& m);
	Matrix operator*(const Matrix& b);
	Matrix operator*(double b);
	Matrix operator/(double b);
	Matrix operator-(const Matrix& m);
	Matrix operator+(const Matrix& m);
	Matrix operator&(const Matrix& b);
	Matrix inv();
	Matrix operator~();
	Matrix operator-();
	Matrix sub(int cs, int cf, int rs, int rf);
	Matrix sub(int s, int f, int inte);
	void setSub(const Matrix& b, int cs, int cf, int rs, int rf);
	void operator*=(double b);
	void operator/=(double b);
	void operator+=(const Matrix& b);
	void operator-=(const Matrix& b);
	double operator()(int row, int col) const;
	void reSet(double data);

	double max();
	double sum();
	double mean();
	double absMean();

	//-getter methods
	int getRow() const;
	int getCol() const;
	void getData(double** data, int& row, int& col) const;
	double* getDataPtr();
	void setData(double* data);               // an existent matrix, change data
	void setData(double* data, int row, int col); //an empty matrix, assign data
	void setData(double data, int row, int col); //change one element of the matrix
	void save(char* file, char* mode);
	int isEmpty();
private:
	void vecLTimesMat(double* ans, double* vec, double* mat, int row, int col);
	void vecRTimesMat(double* ans, double* vec, double* mat, int row, int col);
	void vecTimesVec(double* mat, double* vec);
	void Cramer_Inverse_4x4(double *mat, double *dst);
	int readFile(char* file, int dim, char roc);

	//-private member data
private:
	double* mdp_data;
	int mn_row;
	int mn_col;
	int is_empty;
};

#endif /*MATRIX_H_*/
