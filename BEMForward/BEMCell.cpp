// filename: BEMCell.cpp
// .cpp file used to define class of BEMCell and all its derivative classes

//**************************************************************************************************************
//Inclusions
//**************************************************************************************************************
#include "inclu.h"
#include "external.h"
#include "BEMCell.h"
#include "BEM.h"
//#include "matlab.hpp"

//using namespace std;

//**************************************************************************************************************
//Definitions
//Class:BEMCell
//**************************************************************************************************************
//*******************************************************
//Constructors&destructor
//*******************************************************
BEMCell::BEMCell(int nNodNum, int nDim) :
		m_nNodSeq(1), m_nEleSeq(1), m_nNodNum(nNodNum), m_nDim(nDim)
{
	// pre-collacate storage
	m_dpsinK = new double[m_nNodNum];
	m_dpsinP = new double[m_nNodNum];
	m_dpNodCor = new double[m_nNodNum * m_nDim];
	m_dpPoiCor = new double[m_nDim];
	m_dpNodCon = new double[m_nNodNum * m_nDim];
	m_dpNor = new double[m_nDim];
	m_dpSha = new double[nNodNum];
	// intial value, avoid unpredictable value	
	m_nG = 0;
	for (int j = 0; j < m_nNodNum; j++)
	{
		m_dpsinK[j] = 0;
		m_dpsinP[j] = 0;
		m_dpSha[j] = 0;
	}
	for (int k = 0; k < m_nDim * m_nNodNum; k++)
	{
		m_dpNodCor[k] = 0;
		m_dpNodCon[k] = 0;
	}
	for (int l = 0; l < m_nDim; l++)
	{
		m_dpPoiCor[l] = 0;
		m_dpNor[l] = 0;
	}
}

BEMCell::~BEMCell()
{
	delete[] m_dpsinK;
	delete[] m_dpsinP;
	delete[] m_dpNodCor;
	delete[] m_dpNodCon;
	delete[] m_dpPoiCor;
	delete[] m_dpNor;
	delete[] m_dpSha;
}

//*******************************************************
// function name: setNodSeq
//*******************************************************
void BEMCell::setNodSeq(int nNodSeq, const BEM& bem)
{
	m_nNodSeq = nNodSeq;

	int gpos = (m_nNodSeq - 1) * m_nDim;
	for (int j = 0; j < m_nDim; j++)
		m_dpPoiCor[j] = bem.m_dpNodCor[gpos++];
}

//*******************************************************
// function name: setEleSeq
//*******************************************************
void BEMCell::setEleSeq(int nEleSeq, const BEM& bem)
{
	m_nEleSeq = nEleSeq;

	int seq = 0;
	int pos = 0;
	int gpos = 0;
	for (int i = 0; i < m_nNodNum; i++)
	{
		seq = bem.m_npNodSeq[(m_nEleSeq - 1) * m_nNodNum + i];
		gpos = (seq - 1) * m_nDim;
		for (int j = 0; j < m_nDim; j++)
		{
			m_dpNodCor[pos] = bem.m_dpNodCor[gpos];
//			m_dpNodCon[pos] = bem.m_dpNodCon[gpos];
			pos++;
			gpos++;
		}
	}

	m_nG = 2.0 * bem.m_dpEleAre[nEleSeq - 1];

	pos = (m_nEleSeq - 1) * m_nDim;
	for (int j = 0; j < m_nDim; j++)
		m_dpNor[j] = bem.m_dpEleDer[pos++];

	for (int l = 0; l < m_nNodNum; l++)
	{
		m_dpsinK[l] = 0;
		m_dpsinP[l] = 0;
	}
}

//*******************************************************
// function name: operator =, to avoid = between different derivative classes
//*******************************************************
void BEMCell::operator=(const BEMCell&)
{
}

//**************************************************************************************************************
//Class: Tri3Cell
//**************************************************************************************************************
//*******************************************************
// Constructors/Destructor
//*******************************************************
Tri3Cell::Tri3Cell(int nNodNum, int nDim, int nGauNum, int nsinGauNum) :
		BEMCell(nNodNum, nDim)
{
	genGau(nGauNum, nsinGauNum);

	m_dpsinS = new double[nNodNum];
	for (int i = 0; i < nNodNum; i++)
		m_dpsinS[i] = 0;

}

Tri3Cell::~Tri3Cell()
{
	delete[] m_dpGauPoiCor;
	delete[] m_dpGauPoiWei;
	delete[] m_dpRecGauPoiCor;
	delete[] m_dpRecGauPoiWei;
	delete[] m_dpsinS;
}

//*******************************************************
// function name: genGau
// generate guassian points coordinate and weights according to demanding number of points
// input: nNum --- number of guassian points for triangle
//        nsinNum --number of guassian points for rectangle (singular)
//*******************************************************
void Tri3Cell::genGau(int nNum, int nsinNum)
{
	m_dpGauPoiCor = new double[nNum * 2];
	m_dpGauPoiWei = new double[nNum];
	m_dpRecGauPoiCor = new double[nsinNum];
	m_dpRecGauPoiWei = new double[nsinNum];
	m_nGauNum = nNum;
	m_nSinGauNum = nsinNum;
	switch (nNum)
	{
	case 1:
		m_dpGauPoiCor[0] = m_dpGauPoiCor[1] = 0.333333;
		m_dpGauPoiWei[0] = 1;
		break;
	case 3:
	default:
		m_dpGauPoiCor[0] = 0.166666666666667;
		m_dpGauPoiCor[1] = 0.166666666666667;
		m_dpGauPoiCor[2] = 0.666666666666667;
		m_dpGauPoiCor[3] = 0.166666666666667;
		m_dpGauPoiCor[4] = 0.166666666666667;
		m_dpGauPoiCor[5] = 0.666666666666667;
		/*		m_dpGauPoiCor[0] = 0.5; m_dpGauPoiCor[1] = 0.5;
		 m_dpGauPoiCor[2] = 0.5; m_dpGauPoiCor[3] = 0;
		 m_dpGauPoiCor[4] = 0; m_dpGauPoiCor[5] = 0.5;
		 */
		m_dpGauPoiWei[0] = 0.333333333333333;
		m_dpGauPoiWei[1] = 0.333333333333333;
		m_dpGauPoiWei[2] = 0.333333333333333;
		break;
	case 7:
		m_dpGauPoiCor[0] = 0.101286507323456;
		m_dpGauPoiCor[1] = 0.101286507323456;
		m_dpGauPoiCor[2] = 0.797426985353087;
		m_dpGauPoiCor[3] = 0.101286507323456;
		m_dpGauPoiCor[4] = 0.101286507323456;
		m_dpGauPoiCor[5] = 0.797426985353087;
		m_dpGauPoiCor[6] = 0.470142064105115;
		m_dpGauPoiCor[7] = 0.470142064105115;
		m_dpGauPoiCor[8] = 0.470142064105115;
		m_dpGauPoiCor[9] = 0.059715871789770;
		m_dpGauPoiCor[10] = 0.059715871789770;
		m_dpGauPoiCor[11] = 0.470142064105115;
		m_dpGauPoiCor[12] = 0.333333333333333;
		m_dpGauPoiCor[13] = 0.333333333333333;
		m_dpGauPoiWei[0] = 0.125939180544827;
		m_dpGauPoiWei[1] = 0.125939180544827;
		m_dpGauPoiWei[2] = 0.125939180544827;
		m_dpGauPoiWei[3] = 0.132394152788506;
		m_dpGauPoiWei[4] = 0.132394152788506;
		m_dpGauPoiWei[5] = 0.132394152788506;
		m_dpGauPoiWei[6] = 0.225000000000000;
		break;
	case 13:
		m_dpGauPoiCor[0] = 0.479308067841923;
		m_dpGauPoiCor[1] = 0.260345966079038;
		m_dpGauPoiCor[2] = 0.260345966079038;
		m_dpGauPoiCor[3] = 0.479308067841923;
		m_dpGauPoiCor[4] = 0.2603459660790386;
		m_dpGauPoiCor[5] = 0.260345966079038;
		m_dpGauPoiCor[6] = 0.869739794195568;
		m_dpGauPoiCor[7] = 0.065130102902216;
		m_dpGauPoiCor[8] = 0.065130102902216;
		m_dpGauPoiCor[9] = 0.869739794195568;
		m_dpGauPoiCor[10] = 0.065130102902216;
		m_dpGauPoiCor[11] = 0.065130102902216;
		m_dpGauPoiCor[12] = 0.638444188569809;
		m_dpGauPoiCor[13] = 0.312865496004875;
		m_dpGauPoiCor[14] = 0.638444188569809;
		m_dpGauPoiCor[15] = 0.04869031542532;
		m_dpGauPoiCor[16] = 0.312865496004875;
		m_dpGauPoiCor[17] = 0.638444188569809;
		m_dpGauPoiCor[18] = 0.312865496004875;
		m_dpGauPoiCor[19] = 0.04869031542532;
		m_dpGauPoiCor[20] = 0.04869031542532;
		m_dpGauPoiCor[21] = 0.638444188569809;
		m_dpGauPoiCor[22] = 0.04869031542532;
		m_dpGauPoiCor[23] = 0.312865496004875;
		m_dpGauPoiCor[24] = 0.333333333333333;
		m_dpGauPoiCor[25] = 0.333333333333333;
		m_dpGauPoiWei[12] = -0.149570044467670;
		m_dpGauPoiWei[0] = 0.175615257433204;
		m_dpGauPoiWei[1] = 0.175615257433204;
		m_dpGauPoiWei[2] = 0.175615257433204;
		m_dpGauPoiWei[3] = 0.053347235608839;
		m_dpGauPoiWei[4] = 0.053347235608839;
		m_dpGauPoiWei[5] = 0.053347235608839;
		m_dpGauPoiWei[6] = 0.077113760890257;
		m_dpGauPoiWei[7] = 0.077113760890257;
		m_dpGauPoiWei[8] = 0.077113760890257;
		m_dpGauPoiWei[9] = 0.077113760890257;
		m_dpGauPoiWei[10] = 0.077113760890257;
		m_dpGauPoiWei[11] = 0.077113760890257;
		break;
	}

	switch (nsinNum)
	{
	case 2:
		m_dpRecGauPoiCor[0] = -0.5773502691;
		m_dpRecGauPoiCor[1] = 0.5773502691;
		m_dpRecGauPoiWei[0] = 1;
		m_dpRecGauPoiWei[1] = 1;
		break;
	case 3:
		m_dpRecGauPoiCor[0] = -0.7745966692;
		m_dpRecGauPoiCor[1] = 0;
		m_dpRecGauPoiCor[2] = 0.7745966692;
		m_dpRecGauPoiWei[0] = 0.5555555556;
		m_dpRecGauPoiWei[1] = 0.8888888889;
		m_dpRecGauPoiWei[2] = 0.5555555556;
		break;
	case 4:
		m_dpRecGauPoiCor[0] = -0.3399810436;
		m_dpRecGauPoiCor[1] = -0.8611363116;
		m_dpRecGauPoiCor[2] = 0.3399810436;
		m_dpRecGauPoiCor[3] = 0.8611363116;
		m_dpRecGauPoiWei[0] = 0.6521451549;
		m_dpRecGauPoiWei[1] = 0.3478548451;
		m_dpRecGauPoiWei[2] = 0.6521451549;
		m_dpRecGauPoiWei[3] = 0.3478548451;
		break;
	}
}

//******************************************************
//set the sequence of element, re-initialize
//*******************************************************
void Tri3Cell::setEleSeq(int seq, const BEM& bem)
{
	BEMCell::setEleSeq(seq, bem);
	for (int i = 0; i < m_nNodNum; i++)
		m_dpsinS[i] = 0;
}

//*******************************************************
// Fuction name: norG
// compute |G| which is the result of transform from Cardisan coordinate to natural coordinate
// which is also 2*area of this surface
// 0 x1/ 1 y1 / 2 z1 ; 3 x2/4 y2/5 z2; 6 x3/ 7 y3/ 8 z3; 
//*******************************************************
void Tri3Cell::norG()
{
	double y2y1 = m_dpNodCor[4] - m_dpNodCor[1];
	double y3y1 = m_dpNodCor[7] - m_dpNodCor[1];
	double z2z1 = m_dpNodCor[5] - m_dpNodCor[2];
	double z3z1 = m_dpNodCor[8] - m_dpNodCor[2];
	double x2x1 = m_dpNodCor[3] - m_dpNodCor[0];
	double x3x1 = m_dpNodCor[6] - m_dpNodCor[0];
	m_dpNor[0] = y2y1 * z3z1 - y3y1 * z2z1;   // (y2-y1)(z3-z1) - (y3-y1)(z2-z1)
	m_dpNor[1] = z2z1 * x3x1 - z3z1 * x2x1;    //(z2-z1)(x3-x1) - (z3-z1)(x2-x1)
	m_dpNor[2] = x2x1 * y3y1 - y2y1 * x3x1;    //(x2-x1)(y3-y1) - (y2-y1)(x3-x1)
	m_nG = sqrt(
			m_dpNor[0] * m_dpNor[0] + m_dpNor[1] * m_dpNor[1]
					+ m_dpNor[2] * m_dpNor[2]);
	for (int i = 0; i < m_nDim; i++)
		m_dpNor[i] = m_dpNor[i] / m_nG;
}

//-*******************************************************
// Function name: guaQua
// use Guass Quadrature to compute numeric integral for normal pair
// input: n_shaSeq --- which shape function ( 1, 2, 3 corresponds to N1, N2, and N3)
//        dp_gauCor ---array including gaussian points, for triangle
//        n_guaWeiX, n_guaWeiY --- the weight of guassian points
//        c_flag ---- which coefficient is computed -- H or P
//        type ----- singular or normal
//*******************************************************
//******gauQua for surface integral of u***********************
void Tri3Cell::gauQua_sur(double* data_S, double h, double s, double wr,
		double ws)
{
	m_dpSha[0] = 1 - h - s;
	m_dpSha[1] = h;
	m_dpSha[2] = s;
	for (int i = 0; i < m_nNodNum; i++)
		data_S[i] = data_S[i] + m_dpSha[i] * wr * m_nG * 0.5;

}

//******ordinary**************************************************
void Tri3Cell::gauQua_normal(BEM& bem, int isFirst, double* data_K,
		double* data_P, double h, double s, double wr, double ws)
{
	m_dpSha[0] = 1 - h - s;
	m_dpSha[1] = h;
	m_dpSha[2] = s;

	// x = h1*x1 + h2*x2 + h3*x3;
	double x = m_dpSha[0] * m_dpNodCor[0] + m_dpSha[1] * m_dpNodCor[3]
			+ m_dpSha[2] * m_dpNodCor[6];
	double y = m_dpSha[0] * m_dpNodCor[1] + m_dpSha[1] * m_dpNodCor[4]
			+ m_dpSha[2] * m_dpNodCor[7];
	double z = m_dpSha[0] * m_dpNodCor[2] + m_dpSha[1] * m_dpNodCor[5]
			+ m_dpSha[2] * m_dpNodCor[8];

	//record into the whole matrix, if it is the first sweep
	if (isFirst == 1)
	{
		int pos = bem.m_nGPNum * m_nDim;
		bem.m_dpGps[pos++] = x;
		bem.m_dpGps[pos++] = y;
		bem.m_dpGps[pos] = z;
		bem.m_nGPNum++;
	}

	x = x - m_dpPoiCor[0];
	y = y - m_dpPoiCor[1];
	z = z - m_dpPoiCor[2];

	double r = sqrt(x * x + y * y + z * z);
	double temp_k = -(x * m_dpNor[0] + y * m_dpNor[1] + z * m_dpNor[2]) * wr
			* m_nG * 0.5 / (4 * PI * r * r * r);
	double temp_p = -m_nG * wr * 0.5 / (4 * PI * r);
	// for G type
	for (int i = 0; i < m_nNodNum; i++)
	{
		data_P[i] = data_P[i] + m_dpSha[i] * temp_p;         //seems not correct
		data_K[i] = data_K[i] + m_dpSha[i] * temp_k;
//			data_P[i] = data_P[i] + m_dpSha[i]*wr*m_nG*r*0.5;
	}
}

//****normal for test**********************************
void Tri3Cell::gauQua_normal(double* data_K, double* data_P, double h, double s,
		double wr, double ws)
{
	m_dpSha[0] = 1 - h - s;
	m_dpSha[1] = h;
	m_dpSha[2] = s;

	// x = h1*x1 + h2*x2 + h3*x3;
	double x = m_dpSha[0] * m_dpNodCor[0] + m_dpSha[1] * m_dpNodCor[3]
			+ m_dpSha[2] * m_dpNodCor[6];
	double y = m_dpSha[0] * m_dpNodCor[1] + m_dpSha[1] * m_dpNodCor[4]
			+ m_dpSha[2] * m_dpNodCor[7];
	double z = m_dpSha[0] * m_dpNodCor[2] + m_dpSha[1] * m_dpNodCor[5]
			+ m_dpSha[2] * m_dpNodCor[8];

	x = x - m_dpPoiCor[0];
	y = y - m_dpPoiCor[1];
	z = z - m_dpPoiCor[2];

	double r = sqrt(x * x + y * y + z * z);
	double temp_k = -(x * m_dpNor[0] + y * m_dpNor[1] + z * m_dpNor[2]) * wr
			* m_nG * 0.5 / (4 * PI * r * r * r);
	double temp_p = -m_nG * wr * 0.5 / (4 * PI * r);
	// for G type
	for (int i = 0; i < m_nNodNum; i++)
	{
		data_P[i] = data_P[i] + m_dpSha[i] * temp_p;         //seems not correct
		data_K[i] = data_K[i] + m_dpSha[i] * temp_k;
//			data_P[i] = data_P[i] + m_dpSha[i]*wr*m_nG*r*0.5;
	}
}

//************************singular*************************
void Tri3Cell::gauQua_sin(double* data_K, double* data_P, double h, double s,
		double wr, double ws)
{
	h = 0.5 + 0.5 * h;
	s = 0.5 + 0.5 * s;
	m_dpSha[0] = 1 - h - s * (1 - h);
	m_dpSha[1] = h;
	m_dpSha[2] = s * (1 - h);
	//x = m*x1 + n*x2 + (1-m-n)*x3 / m = a; n = b(1-a)   
	double x = m_dpSha[0] * m_dpNodCor[0] + m_dpSha[1] * m_dpNodCor[3]
			+ m_dpSha[2] * m_dpNodCor[6];
	double y = m_dpSha[0] * m_dpNodCor[1] + m_dpSha[1] * m_dpNodCor[4]
			+ m_dpSha[2] * m_dpNodCor[7];
	double z = m_dpSha[0] * m_dpNodCor[2] + m_dpSha[1] * m_dpNodCor[5]
			+ m_dpSha[2] * m_dpNodCor[8];

	x = x - m_dpPoiCor[0];
	y = y - m_dpPoiCor[1];
	z = z - m_dpPoiCor[2];

	double r = sqrt(x * x + y * y + z * z);
	r = -m_nG * 0.25 * wr * ws * (1 - h) / (4 * PI * r);
//		r = -(x*m_dpNor[0] + y*m_dpNor[1] + z*m_dpNor[2])*wr*ws*m_nG*(1-h)/ (4* PI*r*r*r);
	for (int k = 0; k < m_nNodNum; k++)
	{
		//		data[k] = data[k] + tempt1*dSha[k];                    //seems not correct
		//	   if(k!=0)
		data_P[k] = data_P[k] + r * m_dpSha[k];
	}
}

//*******************************************************
// function name: numInt
// compute numeric integral for each item of H,P and assign its value to corresponding variable
// input: nshaSeq --- which shape function is being dealt with
//*******************************************************
//******pure surface integral of u********************************
void Tri3Cell::surInt()
{
	for (int i = 0; i < m_nGauNum; i++)
		gauQua_sur(m_dpsinS, m_dpGauPoiCor[i * 2], m_dpGauPoiCor[i * 2 + 1],
				m_dpGauPoiWei[i]);
}
//****for the formal BEM computatoin
void Tri3Cell::numInt(BEM& bem, int isFirst)
{
	for (int i = 0; i < m_nGauNum; i++)
		gauQua_normal(bem, isFirst, m_dpsinK, m_dpsinP, m_dpGauPoiCor[i * 2],
				m_dpGauPoiCor[i * 2 + 1], m_dpGauPoiWei[i]);
}

//****for test of single BE integral
void Tri3Cell::numInt()
{
	for (int i = 0; i < m_nGauNum; i++)
		gauQua_normal(m_dpsinK, m_dpsinP, m_dpGauPoiCor[i * 2],
				m_dpGauPoiCor[i * 2 + 1], m_dpGauPoiWei[i]);
}

//*******************************************************
//function name: sinKP 
// construct single vector of PK for each (point, element)pair
// note: need to consider factors such as whether the point is in the element and whether it's related to the shape function
// being processed
//*******************************************************
void Tri3Cell::sinKP(BEM& bem, int isIn, int isFirst)
{
	if (isIn == 0)
		numInt(bem, isFirst);
	else
		sinInt();
}

//*******************************************************
//function name: sigInt
// compute sigular integral for Hij(where i belongs to element, but not = j ),Pii; and Hii(involve Cauchy principle)
//*******************************************************
void Tri3Cell::sinInt()
{
	// change the sequence of nodes within the element such that the sigular points is on the  (1,0) coordinate
	/*	int i = 0;
	 int pos = (m_nEleSeq-1)*m_nNodNum;
	 for ( i=0; i<m_nNodNum; i++)
	 if (bem.m_npNodSeq[pos+i] == m_nNodSeq )
	 break;
	 */
	//check by comparing cor, without the need to access BEM
	int i = 0;
	for (i = 0; i < m_nDim; i++)
	{
		if (m_dpNodCor[3 * i] == m_dpPoiCor[0]
				&& m_dpNodCor[3 * i + 1] == m_dpPoiCor[1]
				&& m_dpNodCor[3 * i + 2] == m_dpPoiCor[2])
			break;
	}

	if (i != 1)
	{
		for (int m = 0; m < m_nDim; m++)
			::swap(m_dpNodCor[m_nDim + m], m_dpNodCor[i * m_nDim + m]);
	}

	//integral for G -- use coordinate trasformation from triangle to rectangle
	//gaussian quadrature
	for (int j = 0; j < m_nSinGauNum; j++)
	{
		for (int k = 0; k < m_nSinGauNum; k++)
		{
			//	gauQua(m_dpsinP,m_dpRecGauPoiCor[j],m_dpRecGauPoiCor[k],m_dpRecGauPoiWei[j],m_dpRecGauPoiWei[k],'G',1);
			//integral for H -- compute Hij with above method,and Hii with "row sum elimination"
			gauQua_sin(m_dpsinK, m_dpsinP, m_dpRecGauPoiCor[j],
					m_dpRecGauPoiCor[k], m_dpRecGauPoiWei[j],
					m_dpRecGauPoiWei[k]);
		}
	}
	//re-change the sequence
//	::swap(m_dpsinP[0],m_dpsinP[pos]);
	if (i != 1)
	{
		::swap(m_dpsinK[1], m_dpsinK[i]);
		::swap(m_dpsinP[1], m_dpsinP[i]);
	}
}

//*******************************************************
//test: whether an integral over a single element is right
//input: cor of 3 apexes
//output: difference between analytic & numeric solutions
//******************************************************
double Tri3Cell::test()
{
	double intNum = 0;
	double intAna = 0;
	int choice;
	cout << "Choose which to test: K (1) or P (2) or Singular (3)" << endl;
	cin >> choice;

	char manual;
	cout << "Manually adjust input coordinates? (y/n)" << endl;
	cin >> manual;
	if (manual == 'y' || manual == 'Y')
	{
		int pos = 0;
		// N1(1,1,0); N2(3,1,0); N3(2,2,0)
		for (int i = 1; i <= m_nNodNum; i++)
		{
			cout << "Node " << i << ": ";
			for (int j = 0; j < m_nDim; j++)
				cin >> m_dpNodCor[pos++];
		}

		cout << "Pls input the coordinate of source point" << endl;
		//(0,0,0) for P, (0,0,1) for K
		for (int j = 0; j < m_nDim; j++)
			cin >> m_dpPoiCor[j];
	}
	else
	{
		if (choice == 2)
		{
			m_dpNodCor[0] = 1;
			m_dpNodCor[1] = 1;
			m_dpNodCor[2] = 0;
			m_dpNodCor[3] = 3;
			m_dpNodCor[4] = 1;
			m_dpNodCor[5] = 0;
			m_dpNodCor[6] = 2;
			m_dpNodCor[7] = 2;
			m_dpNodCor[8] = 0;
			m_dpPoiCor[0] = 0;
			m_dpPoiCor[1] = 0;
			m_dpPoiCor[2] = 0;
		}
		else if (choice == 1)
		{
			m_dpNodCor[0] = 1;
			m_dpNodCor[1] = 1;
			m_dpNodCor[2] = 0;
			m_dpNodCor[3] = 3;
			m_dpNodCor[4] = 1;
			m_dpNodCor[5] = 0;
			m_dpNodCor[6] = 2;
			m_dpNodCor[7] = 2;
			m_dpNodCor[8] = 0;
			m_dpPoiCor[0] = 0;
			m_dpPoiCor[1] = 0;
			m_dpPoiCor[2] = 1;
		}
		else
		{
			m_dpNodCor[0] = 2;
			m_dpNodCor[1] = 1;
			m_dpNodCor[2] = 0;
			m_dpNodCor[3] = 2;
			m_dpNodCor[4] = 2;
			m_dpNodCor[5] = 0;
			m_dpNodCor[6] = 1;
			m_dpNodCor[7] = 2;
			m_dpNodCor[8] = 0;
			m_dpPoiCor[0] = 1;
			m_dpPoiCor[1] = 2;
			m_dpPoiCor[2] = 0;
		}
	}

	norG();
	if (choice == 3)
		sinInt();
	else
		numInt();

	double* valNode = new double[m_nDim];
	double dx, dy, dz;

	if (choice == 1)
	{
		for (int i = 0; i < m_nNodNum; i++)
		{
			dx = m_dpNodCor[m_nDim * i] - m_dpPoiCor[0];
			dy = m_dpNodCor[m_nDim * i + 1] - m_dpPoiCor[1];
			dz = m_dpNodCor[m_nDim * i + 2] - m_dpPoiCor[2];
			valNode[i] = sqrt(dx * dx + dy * dy + dz * dz);
			valNode[i] = valNode[i] * valNode[i] * valNode[i];
			intNum += valNode[i] * m_dpsinK[i];
		}

		intAna = 1 / (4 * PI);

	}
	else if (choice == 2)
	{
		for (int i = 0; i < m_nNodNum; i++)
		{
			dx = m_dpNodCor[m_nDim * i] - m_dpPoiCor[0];
			dy = m_dpNodCor[m_nDim * i + 1] - m_dpPoiCor[1];
			dz = m_dpNodCor[m_nDim * i + 2] - m_dpPoiCor[2];
			valNode[i] = sqrt(dx * dx + dy * dy + dz * dz);
			intNum += valNode[i] * m_dpsinP[i];
		}

		intAna = -DT / (4 * PI * DK);
	}
	else
	{
		for (int i = 0; i < m_nNodNum; i++)
		{
			valNode[i] = 1;
			intNum += valNode[i] * m_dpsinP[i];
		}

		intAna = -DT / (4 * PI * DK) * log(1 + sqrt(2)) * m_nG;
	}

	return (intAna - intNum) / intAna;
}

