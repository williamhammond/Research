#include "inclu.h"
#include "external.h"
#include "Matrix.h"

/**Constructor*/

//************************************************
//ctor(s)/dtor
//**************************************************
//eanble the array 
Matrix::Matrix()
{
   mn_row = 0;
   mn_col = 0;
   mdp_data = NULL;
   is_empty = 1;
}
Matrix::Matrix(double* data, int r, int c):mn_row(r),mn_col(c)
{
	int count = mn_row*mn_col;
	mdp_data = new double[count];
	for(int i=0; i<count; i++)
		mdp_data[i] = data[i];
	is_empty = 0;
}

Matrix::Matrix(int r, int c):mn_row(r),mn_col(c)
{
	int count = mn_row*mn_col;
	mdp_data = new double[count];
	for(int i=0; i<count; i++)
		mdp_data[i] = 0;
	is_empty = 0;
}
//**Copy Constructor*****************
Matrix::Matrix(const Matrix& m)
{
//	m.getData(&mdp_data,mn_row,mn_col);
	mn_row = m.mn_row;
	mn_col = m.mn_col;
	int count = mn_row*mn_col;
	mdp_data = new double[count];
	for(int i=0; i<count; i++)
		mdp_data[i] = m.mdp_data[i];	
	is_empty = 0;
}

//**read from files******************
Matrix::Matrix(char* file, int dim, char roc)
{
	is_empty = readFile(file,dim,roc);
}

//**Deconstructor*******************
Matrix::~Matrix()
{
	if (!is_empty)
		delete []mdp_data;
}


//********************************************************
//access function
//********************************************************
//**get_rows() const*********
//get # of rows of the matrix
//****************************
int Matrix::getRow() const 
{
   return mn_row;
}

//**get_cols() const**********
//get # of columns contained in the matrix.
//***************************
int Matrix::getCol() const
{
   return mn_col;
}

//**getData******************
//write mdp_data into data, also row,cols
//input data is a pointer to the double pointer
//****************************
void Matrix::getData(double** data,int& row, int& col) const
{
	row = this->mn_row;
	col = this->mn_col;
	int count = row*col;
	*data = new double[count];
	for(int i=0; i<count; i++)
		(*data)[i] = this->mdp_data[i];
}


//**getDataPtr****************
//provide access to the pointer to mdp_data
//****************************
double* Matrix::getDataPtr()
{
	return mdp_data;
}
//*******************************
//change the data of the matrix
//*******************************
void Matrix::setData(double* data)
{
	for(int i=0; i<mn_row*mn_col; i++)
		mdp_data[i] = data[i];
}

//*******************************
//assign the data to an empty matrix
//*******************************
void Matrix::setData(double* data, int row, int col)
{
	if (!is_empty)
	{
		cout<<"The matrix is non-empty!!"<<endl;
		::wait();
		exit(0);
	}

	mn_row = row;
	mn_col = col;
	mdp_data = new double[row*col];
	for(int i=0; i<mn_row*mn_col; i++)
		mdp_data[i] = data[i];
	is_empty = 0;
}

//***********************************
//set one element of the matrix
//************************************
void Matrix::setData(double data, int row, int col)
{
	if (is_empty)
	{
		mdp_data = new double[row*col];
		for(int i=0; i<row*col; i++)
			mdp_data[i] = data;
		mn_row = row;
		mn_col = col;
		is_empty = 0;
	}
	else
		mdp_data[col*mn_row+row] = data;
}

//*******************************
//save the data to file
//******************************* 
void Matrix::save(char* file, char* mode)
{
	writeFile(file, mdp_data,sizeof(double),mn_row*mn_col,mode);
//	writeFile(mdp_data,mn_row*mn_col,file,"wb",sizeof(double));
}

//**********************************
//isempty
//************************************
int Matrix::isEmpty()
{
	return is_empty;
}
//*****************************************************
//operator funcitons
//*****************************************************
//operator=: A = B,for copy, not construct (need the constructed object ready)
void Matrix::operator =(const Matrix& m)
{
	if (mn_row != m.getRow() || mn_col != m.getCol())
	{
		cout<<"Error: matrices on lh and rh of '=' must be agree in dimension!"<<endl;
		wait();
		exit(0);
	}
	for(int i=0; i<mn_row*mn_col; i++)
		mdp_data[i] = m.mdp_data[i];
}
//operator(): return A(i,j) with i as row, j as col
double Matrix::operator() (int i, int j) const
{
   return this->mdp_data[j*mn_row+i];
}

//****************************************************
//sub: A(rs,rf,cs,cf) = A(rs:rf,cs:cf);
//****************************************************
Matrix Matrix::sub(int rs,int rf, int cs, int cf)
{
	int row = rf-rs+1;
	int col = cf-cs+1;
	int count = 0, pos = 0;
	Matrix matrix(row,col);
	for(int j=cs; j<=cf; j++)
	{
		pos = j*mn_row + rs;
		for(int i=rs; i<=rf;i++)
			matrix.mdp_data[count++] = mdp_data[pos++];
	}
	return matrix;
}

//****************************************************
//A(s,f,i) = A(s:i:f)
//used for vecotrs
//*****************************************************
Matrix Matrix::sub(int s, int f, int inte )
{
	int count =0; 
	for(int i=s; i<=f; )
	{
		count++;
		i += inte;
	}
	int row,col;
	if (mn_row == 1)
	{
		row = 1;
		col = count;	
	}
	else if (mn_col == 1)
	{
		col = 1;
		row = count;
	}
	else
	{
		cout<<"The matrix is not vecotr!"<<endl;
		wait();
		exit(0);
	}
	Matrix matrix(row,col);
	for(int i = 0; i<row*col; i++)
		matrix.mdp_data[i] = mdp_data[s+i*inte];	
	return matrix;
}

//******************************************************
//setSub: A(rs:rf, cs:cf) = B
//*****************************************************
void Matrix::setSub(const Matrix& b, int rs, int rf, int cs, int cf)
{
	int count = 0, pos = 0;
	for(int j=cs; j<=cf; j++)
	{
		pos = j*mn_row + rs;
		for(int i=rs; i<=rf;i++)
			mdp_data[pos++] = b.mdp_data[count++];
	}
}

//operator+=: A = A + B*******************************
void Matrix::operator+=(const Matrix& m)
{
   if (m.mn_row != mn_row || m.mn_col != mn_col)
   {
	   cout<<"Error: matrix dimension does not agree in addition!"<<endl;
	   wait();
   }
   else
   {
	   for(int i=0; i<mn_row*mn_col;i++)
		   mdp_data[i] += m.mdp_data[i];
   }
}

//operator+=: A = A - B*******************************
void Matrix::operator-=(const Matrix& m)
{
   if (m.mn_row != mn_row || m.mn_col != mn_col)
   {
	   cout<<"Error: matrix dimension does not agree in substraction!"<<endl;
	   wait();
   }
   else
   {
	   for(int i=0; i<mn_row*mn_col;i++)
		   mdp_data[i] -= m.mdp_data[i];
   }

}

//operator*=: A = b*A, b is scalar******************
void Matrix::operator*=(double b)
{
	for(int i=0; i<mn_row*mn_col;i++)
		mdp_data[i] *= b;
}

//operator/=: A = A/b, b is sclar********************
void Matrix::operator/= (double b)
{
	for(int i=0; i<mn_row*mn_col;i++)
		mdp_data[i] /= b;
}

//reSet: A = 0****************************************
void Matrix::reSet(double data)
{
	for(int i=0; i<mn_row*mn_col;i++)
		mdp_data[i] = data;
}

//operator.*: C = A.*B********************************
Matrix Matrix::operator& (const Matrix& m)
{
   if (m.mn_row != mn_row || m.mn_col != mn_col)
   {
	   cout<<"Error: matrix dimension does not agree in dot multiplication!"<<endl;
	   wait();
	   exit(0);
   }
   else
   {
	   Matrix tmp(mn_row,mn_col);
	   for(int i=0; i<mn_row*mn_col;i++)
		   tmp.mdp_data[i] = m.mdp_data[i]*mdp_data[i];
	   
	   return tmp;
   }
}


//operatorP: C = B*b, b is scalar*********************
//*****************************************************
Matrix Matrix::operator *(double b)
{
	Matrix matrix(mn_row,mn_col);
	for (int i=0; i<mn_row*mn_col; i++)
		matrix.mdp_data[i] = b*mdp_data[i];

	return matrix;
}

//operatorP: C = B/b, b is scalar*********************
//*****************************************************
Matrix Matrix::operator /(double b)
{
	Matrix matrix(mn_row,mn_col);
	for (int i=0; i<mn_row*mn_col; i++)
		matrix.mdp_data[i] = mdp_data[i]/b;

	return matrix;
}

//operator*: C = A*B************************************
//consider situations with A is of(1*m) or B is of (m*1). which is much simpler
//*******************************************************
Matrix Matrix::operator* (const Matrix& m)
{
  if(mn_col != m.mn_row)
  {
	  cout<<"Error: matrix dimesion does not agree in multiplication!"<<endl;
	  wait();
	  exit(0);
  }
  else
  {
	  Matrix matrix(mn_row,m.mn_col);
	  // A is 1*c vector
	  if(mn_row == 1)
		  vecLTimesMat(matrix.mdp_data,mdp_data,m.mdp_data,m.mn_row,m.mn_col);
	  else if (m.mn_col == 1)
		  vecRTimesMat(matrix.mdp_data,m.mdp_data,mdp_data,mn_row,mn_col);
	  else
	  {
		  
		   for (int i = 0; i<mn_row; i++)
		   {
			   for (int j=0; j<m.mn_col; j++)
			   {
				   double tmp = 0;
				   for (int k=0; k<mn_col; k++)
				   {
					   tmp += (*this)(i,k)*m(k,j);
					   matrix.mdp_data[j*mn_row+i] = tmp;
				   }
			   }
		   }
	  }
	  return matrix;
  }
}

//operator+: C = A+B************************************
//*******************************************************
Matrix Matrix::operator+ (const Matrix& m)
{
  if(mn_row != m.mn_row || mn_col != m.mn_col)
  {
	  cout<<"Error: matrix dimesion does not agree in addition!"<<endl;
	  wait();
	  exit(0);
  }
  else
  {
	  int count = mn_row*mn_col;
	  Matrix matrix(mn_row,mn_col);
	  for(int i=0; i<count; i++)
		  matrix.mdp_data[i] = mdp_data[i] + m.mdp_data[i];

	  return matrix;
  }
}

//operator+: C = A-B************************************
//*******************************************************
Matrix Matrix::operator- (const Matrix& m)
{
  if(mn_row != m.mn_row || mn_col != m.mn_col)
  {
	  cout<<"Error: matrix dimesion does not agree in substraction!"<<endl;
	  wait();
	  exit(0);
  }
  else
  {
	  int count = mn_row*mn_col;
	  Matrix matrix(mn_row,mn_col);
	  for(int i=0; i<count; i++)
		  matrix.mdp_data[i] = mdp_data[i] - m.mdp_data[i];

	  return matrix;
  }
}


Matrix Matrix::operator- ()
{
	Matrix m(mn_row,mn_col);
	for(int i=0; i<mn_row*mn_col; i++)
		m.mdp_data[i] = -mdp_data[i];
	
	return m;

}



//**** ans = det(A)*********************************
//**************************************************
Matrix Matrix::inv()
{
    Matrix matrix(mn_row,mn_col);
	if (mn_row == 4)
		Cramer_Inverse_4x4(mdp_data,matrix.mdp_data);
/*	else
	{
		engine->addData(mdp_data,"mat",mn_row,mn_col);
		engine->evalString("mat = inv(mat)");
		engine->getData(data,"mat");
	}
*/	
	return matrix;
}

//operator~: A = transpose(A)*******************************
//******************************************************
Matrix Matrix::operator~ ()
{
	Matrix matrix(mn_col,mn_row);
	for(int i=0; i<mn_col; i++)
		for(int j=0; j<mn_row; j++)
			matrix.mdp_data[j*mn_col+i] = (*this)(j,i);
   return matrix;
}

//max:return the maximum element
//*******************************************************
double Matrix::max()
{
	double tmp = mdp_data[0];
	for(int i=0; i<mn_row*mn_col; i++)
	{
		if(tmp<mdp_data[i])
			tmp = mdp_data[i];
	}
	return tmp;
}

//sum:return the sum of all elements
//*******************************************************
double Matrix::sum()
{
	double tmp = 0;
	for(int i=0; i<mn_row*mn_col; i++)
		tmp += mdp_data[i];
	return tmp;
}

//mean:return mean of all elements
//*******************************************************
double Matrix::mean()
{
	double tmp = 0;
	for(int i=0; i<mn_row*mn_col; i++)
		tmp += mdp_data[i];
	return tmp/(mn_row*mn_col);
}

//mean:return mean of all abs(elements)
//*******************************************************
double Matrix::absMean()
{
	double tmp = 0;
	for(int i=0; i<mn_row*mn_col; i++)
		tmp += abs(mdp_data[i]);
	return tmp/(mn_row*mn_col);
}
//********************************************************
//ans = vec * mat, where ans, vec is vecs of 1*row, mat is row*col
//*******************************************************
void Matrix::vecLTimesMat(double* ans,double* vec,double* mat, int row, int col)
{
	::reSet(ans,col);
	int count = 0;
	for(int i=0; i<col;i++)
		for(int j = 0; j<row; j++)
			ans[i] += vec[j]*mat[count++];
 
}

//********************************************************
//ans = mat*vec, where ans, vec is vecs of col*1, mat is row*col
//*******************************************************
void Matrix::vecRTimesMat(double* ans,double* vec,double* mat, int row, int col)
{
	::reSet(ans,row);
	int count = 0;
	for(int i=0; i<col;i++)
		for(int j = 0; j<row; j++)
			ans[j] += vec[i]*mat[count++];
 
}


void Matrix::Cramer_Inverse_4x4(double *mat, double *dst)
{
	double tmp[12]; /* temp array for pairs */
	double src[16]; /* array of transpose source matrix */
	double det; /* determinant */
	/* transpose matrix */
	for (int i = 0; i < 4; i++) {
	src[i] = mat[i*4];
	src[i + 4] = mat[i*4 + 1];
	src[i + 8] = mat[i*4 + 2];
	src[i + 12] = mat[i*4 + 3];
	}
	/* calculate pairs for first 8 elements (cofactors) */
	tmp[0] = src[10] * src[15];
	tmp[1] = src[11] * src[14];
	tmp[2] = src[9] * src[15];
	tmp[3] = src[11] * src[13];
	tmp[4] = src[9] * src[14];
	tmp[5] = src[10] * src[13];
	tmp[6] = src[8] * src[15];
	tmp[7] = src[11] * src[12];
	tmp[8] = src[8] * src[14];
	tmp[9] = src[10] * src[12];
	tmp[10] = src[8] * src[13];
	tmp[11] = src[9] * src[12];
	/* calculate first 8 elements (cofactors) */
	dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
	dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
	dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
	dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
	dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
	dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
	dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
	dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
	dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
	dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
	dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
	dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
	dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
	dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
	dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
	dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
	/* calculate pairs for second 8 elements (cofactors) */
	tmp[0] = src[2]*src[7];
	tmp[1] = src[3]*src[6];
	tmp[2] = src[1]*src[7];
	tmp[3] = src[3]*src[5];
	tmp[4] = src[1]*src[6];
	tmp[5] = src[2]*src[5];
	tmp[6] = src[0]*src[7];
	tmp[7] = src[3]*src[4];
	tmp[8] = src[0]*src[6];
	tmp[9] = src[2]*src[4];
	tmp[10] = src[0]*src[5];
	tmp[11] = src[1]*src[4];
	/* calculate second 8 elements (cofactors) */
	dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
	dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
	dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
	dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
	dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
	dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
	dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
	dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
	dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
	dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
	dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
	dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
	dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
	dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
	dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
	dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
	/* calculate determinant */
	det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
	if(det<DET_INVERSE)
	{
		cout<<"Error: support nodes not enought!"<<endl;
		wait();
	}
	/* calculate matrix inverse */
	det = 1/det;
//	cout<<"det = "<<det<<endl;
	for (int j = 0; j < 16; j++)
		dst[j] *= det;
}

//*************************************************
//read in files, if file exist, construct matrix, return 0
//otherwise, cnstruct empty matrix,return 1
//*************************************************
int  Matrix::readFile(char* file, int dim, char roc)
{
	FILE* hfil = fopen(file,"rb");
	if(hfil==NULL)
	{
		mn_row = 0;
		mn_col = 0;
		mdp_data = NULL;
		return 1;
	}
	else
	{
		if (roc == 'r')
			mn_row = dim;
		else
			mn_col = dim;

		fseek(hfil,0,SEEK_END);
        int filSiz = ftell(hfil);
		rewind(hfil);
        //get the size(rowDim*colDim)of input matrix
		int datByt = filSiz/sizeof(double);
		if (roc == 'r')
			mn_col = datByt/dim;
		else
			mn_row = datByt/dim;
		mdp_data = new double[datByt];
		fread(mdp_data,sizeof(double),datByt,hfil); // read out colDim data for every row, each data size is specified by"size"
		fclose(hfil);
		return 0;
	}
}

//**********************************************
//two matrix A = A+scale*B
//*************************************************
/*
void matAdd(double* a, double* b, double scale,  int len)
{
	for(int i=0; i<len; i++)
		a[i] += (scale*b[i]);
}
*/
//**************************************************
//matrix B, vector b, augment B to [B,scale*b]
//**************************************************
/*
void matAug(double* B, double* b, double scale, int len)
{
	for(int i=0; i<len; i++)
         B[i] = scale*b[i];
}
*/



