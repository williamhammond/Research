// Filename :BEM.cpp
// definition of BEM class & external functions

//**************************************************************************************************************
//inclusions
//**************************************************************************************************************
#include"inclu.h"
#include"external.h"
#include"BEMCell.h"
#include"BEM.h"
#include"MFree.h"
//#include "matlab.hpp"
//**************************************************************************************************************
// Class : BEM
// Function: 1, generate the whole BEM mesh ( input interface: read in the file)
//           2, process the whole method ( inner, call BEMClass object)
//           3, give out the final result, including solution, K,P matrix, transfer matrix for inverse
//**************************************************************************************************************
//*******************************************************
//ctor/dtor
//*******************************************************
BEM::BEM(char** filNamBEM, char** filNamMFree, CellType celTyp)
{
	// generate the BEMClass object according to required CelTyp
	switch (celTyp)
	{
	case Tri3:
		m_ptrBEMEle = new Tri3Cell(3, 3, GAUPOI, SINGAUPOI); // can be changed later for flexible control
		break;
	case Lin2:
		break;
	case Rec4:
		break;
	case Tri6:
		break;
		//....
	}

	// readfile to generate the mesh & allocate storage
	// read in data for surface nodes
	m_dpNodCor = NULL;
	m_npNodSeq = NULL;
	m_npForbidden = NULL;
	readFile(filNamBEM[0], &m_dpNodCor, sizeof(double), m_ptrBEMEle->m_nDim,
			m_nNodNum);
	readFile(filNamBEM[1], &m_npNodSeq, sizeof(int), m_ptrBEMEle->m_nNodNum,
			m_nEleNum);
	if (FORBIDDEN)
		readFile(filNamBEM[2], &m_npForbidden, sizeof(int), 1, m_nForNum);
	else
	{
		m_nForNum = 0;
		m_npForbidden = NULL;
	}
//	readFil(filNamBEM[2],&m_dpNodCon,sizeof(double),m_ptrBEMEle->m_nDim,m_nNodNum);

	init(filNamMFree);
}
//****no MFree part, i.e.Laplacian equation*********************************************************
BEM::BEM()
{
}

BEM::~BEM()
{
	//delete object pointer
	delete m_ptrBEMEle;
	delete m_ptrMFree;
	//delete 2-D matrices
	delete[] m_npNodSeq;
	delete[] m_dpNodCor;
//	delete []m_dpNodCon;
	delete[] m_dpEleDer;
	delete[] m_dpEleAre;
	delete[] m_dpRefCor;
	delete[] m_dpK;
	delete[] m_dpP;
	delete[] m_dpK_arranged;
	delete[] m_dpP_arranged;
	delete[] m_dpGps;
	if (m_npForbidden)
		delete[] m_npForbidden;
	if (m_npDiscarded)
		delete[] m_npDiscarded;
}

void BEM::init(char** filNamMFree)
{
	int count = 0;
	m_ptrMFree = new MFree(BG_INTERVAL, filNamMFree, *this);
//	m_dpRefCor = new double[m_ptrBEMEle->m_nDim];
//	for(int j=0; j<m_ptrBEMEle->m_nDim; j++)
//		m_dpRefCor[j] = 0;
	//preallocate storage
	//for BEM
	count = m_nEleNum * m_ptrBEMEle->m_nDim;
	m_dpEleDer = new double[count];
	for (int i = 0; i < count; i++)
		m_dpEleDer[i] = 0;
	m_dpEleAre = new double[m_nEleNum];
	for (int k = 0; k < m_nEleNum; k++)
		m_dpEleAre[k] = 0;
	varCor.x = varCor.y = varCor.z = 0;

	count = m_nNodNum * m_nNodNum;
	m_dpP = new double[count];
	m_dpK = new double[count];
	for (int l = 0; l < count; l++)
	{
		m_dpP[l] = 0;
		m_dpK[l] = 0;
	}

	cout << "Allocate for P, K arraged" << endl;

	count = (m_nNodNum - m_nForNum) * (m_nNodNum - m_nForNum);
	m_dpP_arranged = new double[count];
	m_dpK_arranged = new double[count];
	for (int h = 0; h < count; h++)
	{
		m_dpP_arranged[h] = 0;
		m_dpK_arranged[h] = 0;
	}

	if (m_nForNum)
	{
		m_nDiscarded = 2 * m_nNodNum * m_nForNum; //# elements to be discarded in the matrice
		//including repetition, didn't concern it for 
		//the sake of simplicity
		m_npDiscarded = new int[m_nDiscarded];
		count = 0;
		int forSeq;
		for (int a = 0; a < m_nForNum; a++)
		{
			forSeq = m_npForbidden[a] - 1;
			//the k-th row
			for (int b = 0; b < m_nNodNum; b++)
				m_npDiscarded[count++] = forSeq * m_nNodNum + b;
			//the k-th col
			for (int c = 0; c < m_nNodNum; c++)
				m_npDiscarded[count++] = c * m_nNodNum + forSeq;

		}
		if (count != m_nDiscarded)
		{
			cout << "Error in determining discarded elements!" << endl;
			exit(0);
		}
	}
	else
	{
		m_nDiscarded = 0;
		m_npDiscarded = NULL;
	}

	//only record in the first sweep, since each repetion have the same gps except differneces
	//in singular elements. (No need to record down sigular ones)
	count = m_nEleNum * GAUPOI * m_ptrBEMEle->m_nDim;
	m_dpGps = new double[count];
	for (int m = 0; m < count; m++)
		m_dpGps[m] = 0;
	m_nGPNum = 0;

	m_bisFirst = true;

}

//*******************************************************
//function name: isBel
//judge whether node i is now belonging to element j
//*******************************************************
int BEM::isBel(int nNodeSeq, int nEleSeq)
{
	int pos = (nEleSeq - 1) * m_ptrBEMEle->m_nNodNum;
	for (int i = 0; i < m_ptrBEMEle->m_nNodNum; i++)
	{
		if (nNodeSeq == m_npNodSeq[pos++])
			return 1;
	}
	return 0;
}

//*******************************************************
//isForbidden
//judge whether the node (nodSeq) is used in BEM
//******************************************************
bool BEM::isForbidden(int nodSeq)
{
	if (FORBIDDEN)
	{
		for (int i = 0; i < m_nForNum; i++)
		{
			if (nodSeq == m_npForbidden[i])
				return true;
		}
		return false;
	}
	else
		return false;
}

//*************************************************************
//on ordinary nodes (not the first loop)
//************************************************************
void BEM::intPK(int nNodSeq, int neleSeq)
{
	//generate single K,P
	m_ptrBEMEle->sinKP(*this, isBel(nNodSeq, neleSeq), 0);
//	cout<<m_ptrBEMEle->m_dpsinK[0]<<"  "<<m_ptrBEMEle->m_dpsinK[1]<<"  "<<m_ptrBEMEle->m_dpsinK[2]<<endl;
	// assemble into integral K,P
	int pos = (neleSeq - 1) * m_ptrBEMEle->m_nNodNum;
	int col = (nNodSeq - 1) * m_nNodNum;
	int seq = 0;
//	wriFil(".\\heart836torso700\\output\\eleK.bin",m_ptrBEMEle->m_dpsinK, sizeof(double), 3, "a+b");
	for (int i = 0; i < m_ptrBEMEle->m_nNodNum; i++)
	{
		seq = m_npNodSeq[pos++];
//        wriFil(".\\heart836torso700\\output\\ele.bin",&seq, sizeof(int), 1, "a+b");
		m_dpK[col + seq - 1] += m_ptrBEMEle->m_dpsinK[i];
		m_dpP[col + seq - 1] += m_ptrBEMEle->m_dpsinP[i];
	}
//	double* ptr = mxGetPr(m_dpK.GetData());
//	wriFil(".\\heart836torso700\\test\\testK.bin",ptr,sizeof(double),(m_nNodNum+1)*m_nNodNum,"wb");
}

//********************************************************
//rearrangePK
//get rid of forbidden nodes
//*********************************************************
void BEM::rearrangePK()
{
	int count = m_nNodNum * m_nNodNum;
	int pos = 0;
	cout << endl;
	cout << "************ Rearranging **********" << endl;
	if (FORBIDDEN)
	{
		for (int i = 0; i < count; i++)
		{
//        cout<<i<<" "<<endl;
			if (!isDiscarded(i))
			{
				m_dpK_arranged[pos] = m_dpK[i];
				m_dpP_arranged[pos++] = m_dpP[i];
			}
		}
		count = m_nNodNum - m_nForNum;
		if (pos != count * count)
		{
			cout << "Errors during the rearrangement of P&K !" << endl;
			exit(0);
		}
	}
	else
	{
		for (int i = 0; i < count; i++)
		{
			m_dpK_arranged[pos] = m_dpK[i];
			m_dpP_arranged[pos++] = m_dpP[i];
		}
	}
}

//******************************************************************
//isDiscarded
//whether i is in the k_th row or col of P,K matrices, where k corresponds
//to the seq of forbidden nodes.
//if so, this is removed from the final P,Kmatrices
//*******************************************************************
bool BEM::isDiscarded(int seq)
{
	for (int i = 0; i < m_nDiscarded; i++)
	{
		if (seq == m_npDiscarded[i])
			return true;
	}
	return false;
}

//*******************************************************************************
//visit function
//******************************************************************************
double* BEM::getBEMK()
{
	return m_dpK_arranged;
}

double* BEM::getBEMP()
{
	return m_dpP_arranged;
}

int BEM::getNumSur()
{
	return m_nNodNum - m_nForNum;
}

int BEM::getNumVol()
{
	return m_ptrMFree->m_nNodNum;
}
//************************************************************************************************************************
//BEM_source (for Poisson)
//************************************************************************************************************************
//ctor/dtor
BEM_source::BEM_source(char** filNamBEM, char**filNamMFree, CellType type) :
		BEM(filNamBEM, filNamMFree, type)
{
	m_dpRefCor = new double[m_ptrBEMEle->m_nDim];
	for (int j = 0; j < m_ptrBEMEle->m_nDim; j++)
		m_dpRefCor[j] = 0;
	//tetra for torso************
	m_dpTetCor = NULL;
	m_npTet = NULL;
	readFile(filNamBEM[4], &m_dpTetCor, sizeof(double), m_ptrBEMEle->m_nDim,
			m_nTetNodNum);
	readFile(filNamBEM[5], &m_npTet, sizeof(int), 4, m_nTetNum);

	m_dpKa = new double[m_nNodNum];
	for (int l = 0; l < m_nNodNum; l++)
		m_dpKa[l] = 0;
	m_dpKa_arranged = new double[m_nNodNum - m_nForNum];
	for (int k = 0; k < m_nNodNum - m_nForNum; k++)
		m_dpKa_arranged[k] = 0;
	for (int m = 0; m < m_nNodNum; m++)
	{
		for (int n = 0; n < m_ptrBEMEle->m_nDim; n++)
			m_dpRefCor[n] += m_dpNodCor[m_ptrBEMEle->m_nDim * m + n];
	}
	for (int q = 0; q < m_ptrBEMEle->m_nDim; q++)
		m_dpRefCor[q] /= m_nNodNum;

	//compute the center of each element, used to visualize "surface normal field"
	int count = m_ptrBEMEle->m_nDim * m_nEleNum;
	double* cor_center = new double[count];
	int node1, node2, node3, pos;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	for (int i = 0; i < m_nEleNum; i++)
	{
		pos = i * m_ptrBEMEle->m_nNodNum;
		node1 = m_npNodSeq[pos++] - 1;
		node2 = m_npNodSeq[pos++] - 1;
		node3 = m_npNodSeq[pos++] - 1;
		pos = node1 * m_ptrBEMEle->m_nDim;
		x1 = m_dpNodCor[pos++];
		y1 = m_dpNodCor[pos++];
		z1 = m_dpNodCor[pos++];
		pos = node2 * m_ptrBEMEle->m_nDim;
		x2 = m_dpNodCor[pos++];
		y2 = m_dpNodCor[pos++];
		z2 = m_dpNodCor[pos++];
		pos = node3 * m_ptrBEMEle->m_nDim;
		x3 = m_dpNodCor[pos++];
		y3 = m_dpNodCor[pos++];
		z3 = m_dpNodCor[pos++];
		pos = i * m_ptrBEMEle->m_nDim;
		cor_center[pos++] = (x1 + x2 + x3) / 3.0;
		cor_center[pos++] = (y1 + y2 + y3) / 3.0;
		cor_center[pos++] = (z1 + z2 + z3) / 3.0;
	}

	writeFile(filNamBEM[3], cor_center, sizeof(double), count, "wb");
}

BEM_source::~BEM_source()
{
	delete[] m_dpTetCor;
	delete[] m_npTet;
	delete[] m_dpKa;
	delete[] m_dpKa_arranged;
}

//*******************************************************
//isOutward: test the firstly computed Nor is outward or inward
//if element on the torso surface, easy, just use a ref point within torso
//if element on heart suface, use norm to generate a point, then use tet to test whether that point is in 
//*******************************************************
int BEM_source::isOutward()
{
	if (OUT_BOUNDARY == 'T')
	{
		if (!IS_CONVEX)
		{
			varCor.x = m_ptrBEMEle->m_dpNodCor[0]
					- SCALE_T * m_ptrBEMEle->m_dpNor[0];
			varCor.y = m_ptrBEMEle->m_dpNodCor[1]
					- SCALE_T * m_ptrBEMEle->m_dpNor[1];
			varCor.z = m_ptrBEMEle->m_dpNodCor[2]
					- SCALE_T * m_ptrBEMEle->m_dpNor[2];
			if (isoutTetra(varCor, m_dpTetCor, m_npTet, m_nTetNum))
				return 0;
			else
				return 1;
		}
		else
		//not always right, especially in neck, arm region
		{
			double dx = m_ptrBEMEle->m_dpNodCor[0] - m_dpRefCor[0];
			double dy = m_ptrBEMEle->m_dpNodCor[1] - m_dpRefCor[1];
			double dz = m_ptrBEMEle->m_dpNodCor[2] - m_dpRefCor[2];
			if (dx * m_ptrBEMEle->m_dpNor[0] + dy * m_ptrBEMEle->m_dpNor[1]
					+ dz * m_ptrBEMEle->m_dpNor[2] > 0)
				return 1;
			else
				return 0;
		}
	}
	else
	{
		varCor.x = m_ptrBEMEle->m_dpNodCor[0]
				- SCALE_T * m_ptrBEMEle->m_dpNor[0];
		varCor.y = m_ptrBEMEle->m_dpNodCor[1]
				- SCALE_T * m_ptrBEMEle->m_dpNor[1];
		varCor.z = m_ptrBEMEle->m_dpNodCor[2]
				- SCALE_T * m_ptrBEMEle->m_dpNor[2];
		if (m_ptrMFree->isout(varCor))
			return 0;
		else
			return 1;
	}

}

//*******************************************************
// eleNor
// record normal derivative & area for every triangular surface of BEM
//*******************************************************
void BEM_source::eleNor(int neleSeq)
{
	m_ptrBEMEle->norG();
	if (!isOutward())
	{
		for (int j = 0; j < m_ptrBEMEle->m_nDim; j++)
			m_ptrBEMEle->m_dpNor[j] = -m_ptrBEMEle->m_dpNor[j];
	}
	int pos = (neleSeq - 1) * m_ptrBEMEle->m_nDim;
	for (int i = 0; i < m_ptrBEMEle->m_nDim; i++)
		m_dpEleDer[pos++] = m_ptrBEMEle->m_dpNor[i];
	m_dpEleAre[neleSeq - 1] = absolute(m_ptrBEMEle->m_nG) / 2.0;
}

//*************************************************************
//intPK: for the first loop, node 1
//*************************************************************
void BEM_source::intPK_first(int nnodSeq, int neleSeq)
{
	eleNor(neleSeq);
	//generate single K,P
	m_ptrBEMEle->sinKP(*this, isBel(nnodSeq, neleSeq), 1);

	m_ptrBEMEle->surInt();
//	cout<<m_ptrBEMEle->m_dpsinK[0]<<"  "<<m_ptrBEMEle->m_dpsinK[1]<<"  "<<m_ptrBEMEle->m_dpsinK[2]<<endl;
	// assemble into integral K,P
	int pos = (neleSeq - 1) * m_ptrBEMEle->m_nNodNum;
	int seq = 0;
//	wriFil(".\\heart836torso700\\output\\eleK.bin",m_ptrBEMEle->m_dpsinK, sizeof(double), 3, "a+b");
	for (int i = 0; i < m_ptrBEMEle->m_nNodNum; i++)
	{
		seq = m_npNodSeq[pos++];
//       wriFil(".\\heart836torso700\\output\\ele.bin",&seq, sizeof(int), 1, "a+b");
		m_dpK[(nnodSeq - 1) * m_nNodNum + seq - 1] += m_ptrBEMEle->m_dpsinK[i];
		m_dpP[(nnodSeq - 1) * m_nNodNum + seq - 1] += m_ptrBEMEle->m_dpsinP[i];
		m_dpKa[seq - 1] += m_ptrBEMEle->m_dpsinS[i]; // m_dpEleAre[neleSeq-1];
	}

//	double* ptr = mxGetPr(m_dpK.GetData());
//	wriFil(".\\heart836torso700\\test\\testK.bin",ptr,sizeof(double),(m_nNodNum+1)*m_nNodNum,"wb");
}

//*******************************************************
//BEMCompute
//final interface with external application, generate final P,K,B,A and Transfer matrix
//*******************************************************
void BEM_source::BEMCompute(char** filOutput)
{
	clock_t s, f;
	double d;
	s = clock();
	//loop for each surface node
	cout << "********BEM*****************************************" << endl;
	for (int i = 1; i <= m_nNodNum; i++)
	{
		//loop for each element in surface
		if (!isForbidden(i))
		{
			cerr << i;
			m_ptrBEMEle->setNodSeq(i, *this);
			if (m_bisFirst)
			{
				for (int j = 1; j <= m_nEleNum; j++)
				{
					m_ptrBEMEle->setEleSeq(j, *this);
					intPK_first(i, j);
				}
				m_bisFirst = false;
			}
			else
			{
				for (int j = 1; j <= m_nEleNum; j++)
				{
					m_ptrBEMEle->setEleSeq(j, *this);
					intPK(i, j);
				}
			}
			//		test("K_step_source.bin",&m_dpK[(i-1)*m_nNodNum],sizeof(double),m_nNodNum,"wb");
		}
	}

	writeFile(
			"/Users/maomaowlw/Research/Mywork/Data/Geometry/CESC11/Healthy/Processed/Ref/heart_vs.normal",
			m_dpEleDer, sizeof(double), m_nEleNum * m_ptrBEMEle->m_nDim, "wb");
	/* 	test("BEM_outward_normal.bin", m_dpEleDer,sizeof(double), m_nEleNum*m_ptrBEMEle->m_nDim, "wb");
	 test("BEM_gp.cor",m_dpGps,sizeof(double),m_nGPNum*m_ptrBEMEle->m_nDim,"wb");
	 test("original_BEMK.bin",m_dpK,sizeof(double),m_nNodNum*m_nNodNum,"wb");
	 test("original_BEMP.bin",m_dpP,sizeof(double),m_nNodNum*m_nNodNum,"wb");
	 test("original_Ka.bin",m_dpKa,sizeof(double),m_nNodNum,"wb");
	 */
	rearrangePK();
	// get Hii from row sum elimination
	int k, pos = 0;
	int num = m_nNodNum - m_nForNum;
	double temp = 0;
	for (k = 0; k < num; k++)
	{
		temp = 0;
		pos = k * num;
		for (int l = 0; l < num; l++)
		{
			if (l != k)
				temp += m_dpK_arranged[pos + l];
		}
		m_dpK_arranged[pos + k] = -temp;
	}

	//rescale the area
	temp = 0;
	for (int l = 0; l < m_nEleNum; l++)
		temp += m_dpEleAre[l];
	for (k = 0; k < num; k++)
		m_dpKa_arranged[k] /= temp;
	//write out
	writeFile(filOutput[0], m_dpK_arranged, sizeof(double), num * num, "wb");
	writeFile(filOutput[1], m_dpP_arranged, sizeof(double), num * num, "wb");
	writeFile(filOutput[2], m_dpKa_arranged, sizeof(double), num, "wb");
	f = clock();
	d = double(f - s) / CLOCKS_PER_SEC;
	::cout << "time for BEM" << d << endl;
	//MFree compute B
	if (!FLAG)
		m_ptrMFree->assTrans(*this, filOutput[3]);
//	wriFile(filOutput[3],m_ptrMFree->Trans,"wb");
	/*	else
	 {
	 double* con = NULL;
	 readFil(".\\cube64\\Parameter_normal.bin",&con,sizeof(double),2,m_ptrMFree->m_nNodNum);
	 m_ptrMFree->assTrans(con,*this);
	 }*/
//	cout<<"partially finishing"<<endl;
}

//***************************************************
//rearrange
//***************************************************
void BEM_source::rearrangePK()
{
	BEM::rearrangePK();
	int pos = 0;
	if (FORBIDDEN)
	{
		for (int j = 0; j < m_nNodNum; j++)
		{
			if (!isForbidden(j + 1))
				m_dpKa_arranged[pos++] = m_dpKa[j];
		}
		if (pos != m_nNodNum - m_nForNum)
		{
			cout << "Errors during the rearrangement of Ka !" << endl;
			exit(0);
		}
	}
	else
	{
		for (int j = 0; j < m_nNodNum; j++)
			m_dpKa_arranged[pos++] = m_dpKa[j];
	}
}
//***************************************************
//additional visit function
//***************************************************		
double* BEM_source::getBEMKa()
{
	return m_dpKa_arranged;
}

double* BEM_source::getMFree()
{
	return m_ptrMFree->m_dpTrans_arranged;
}

//*********************************************************************************************************************
//BEM no source, pure BEM
//********************************************************************************************************************
//ctor/dtor
BEM_nosource::BEM_nosource(char** filNamBEM, char** filNamMFree, CellType type)
{
	switch (type)
	{
	case Tri3:
		m_ptrBEMEle = new Tri3Cell(3, 3, GAUPOI, SINGAUPOI); // can be changed later for flexible control
		break;
		//....
	}
	double* cor_inner = NULL, *cor_outer = NULL;
	int* seq_inner = NULL, *seq_outer = NULL, *for_inner = NULL, *for_outer =
			NULL;
	readFile(filNamBEM[0], &cor_inner, sizeof(double), m_ptrBEMEle->m_nDim,
			m_nNodNumInner);
	readFile(filNamBEM[1], &seq_inner, sizeof(int), m_ptrBEMEle->m_nNodNum,
			m_nEleNumInner);
	if (FORBIDDEN)
		readFile(filNamBEM[2], &for_inner, sizeof(int), 1, m_nForNumInner);
	else
		m_nForNumInner = 0;
	readFile(filNamBEM[4], &cor_outer, sizeof(double), m_ptrBEMEle->m_nDim,
			m_nNodNumOuter);
	readFile(filNamBEM[5], &seq_outer, sizeof(int), m_ptrBEMEle->m_nNodNum,
			m_nEleNumOuter);
	if (FORBIDDEN)
		readFile(filNamBEM[6], &for_outer, sizeof(int), 1, m_nForNumOuter);
	else
		m_nForNumOuter = 0;

	m_nNodNum = m_nNodNumInner + m_nNodNumOuter;
	m_dpNodCor = new double[m_ptrBEMEle->m_nDim * m_nNodNum];
	int count = m_ptrBEMEle->m_nDim * m_nNodNumInner;
	for (int h = 0; h < count; h++)
		m_dpNodCor[h] = cor_inner[h];
	for (int j = count; j < m_ptrBEMEle->m_nDim * m_nNodNum; j++)
		m_dpNodCor[j] = cor_outer[j - count];
	delete[] cor_inner;
	delete[] cor_outer;

	m_nEleNum = m_nEleNumInner + m_nEleNumOuter;
	m_npNodSeq = new int[m_ptrBEMEle->m_nNodNum * m_nEleNum];
	count = m_ptrBEMEle->m_nNodNum * m_nEleNumInner;
	for (int l = 0; l < count; l++)
		m_npNodSeq[l] = seq_inner[l];
	for (int k = count; k < m_ptrBEMEle->m_nNodNum * m_nEleNum; k++)
		m_npNodSeq[k] = m_nNodNumInner + seq_outer[k - count];
	delete[] seq_inner;
	delete[] seq_outer;
	m_nForNum = m_nForNumInner + m_nForNumInner;
	if (m_nForNum)
	{
		m_npForbidden = new int[m_nForNum];
		for (int t = 0; t < m_nForNumInner; t++)
			m_npForbidden[t] = for_inner[t];
		for (int r = m_nForNumInner; r < m_nForNum; r++)
			m_npForbidden[r] = for_outer[r - m_nForNumInner] + m_nNodNumInner;
		if (for_inner)
			delete[] for_inner;
		if (for_outer)
			delete[] for_outer;
	}
	else
		m_npForbidden = NULL;

	//write out, to facilitate later visulization.
	writeFile(filNamBEM[7], m_dpNodCor, sizeof(double),
			m_ptrBEMEle->m_nDim * m_nNodNum, "wb");
	writeFile(filNamBEM[8], m_npNodSeq, sizeof(int),
			m_ptrBEMEle->m_nNodNum * m_nEleNum, "wb");
	if (m_npForbidden)
		writeFile(filNamBEM[9], m_npForbidden, sizeof(int), m_nForNum, "wb");

	//compute the center of each element, used to visualize "surface normal field"
	count = m_ptrBEMEle->m_nDim * m_nEleNum;
	double* cor_center = new double[count];
	int node1, node2, node3, pos;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	for (int i = 0; i < m_nEleNum; i++)
	{
		pos = i * m_ptrBEMEle->m_nNodNum;
		node1 = m_npNodSeq[pos++] - 1;
		node2 = m_npNodSeq[pos++] - 1;
		node3 = m_npNodSeq[pos++] - 1;
		pos = node1 * m_ptrBEMEle->m_nDim;
		x1 = m_dpNodCor[pos++];
		y1 = m_dpNodCor[pos++];
		z1 = m_dpNodCor[pos++];
		pos = node2 * m_ptrBEMEle->m_nDim;
		x2 = m_dpNodCor[pos++];
		y2 = m_dpNodCor[pos++];
		z2 = m_dpNodCor[pos++];
		pos = node3 * m_ptrBEMEle->m_nDim;
		x3 = m_dpNodCor[pos++];
		y3 = m_dpNodCor[pos++];
		z3 = m_dpNodCor[pos++];
		pos = i * m_ptrBEMEle->m_nDim;
		cor_center[pos++] = (x1 + x2 + x3) / 3.0;
		cor_center[pos++] = (y1 + y2 + y3) / 3.0;
		cor_center[pos++] = (z1 + z2 + z3) / 3.0;
	}

	writeFile(filNamBEM[10], cor_center, sizeof(double), count, "wb");
	delete[] cor_center;

	init(filNamMFree);

	m_dpRefCor = new double[m_ptrBEMEle->m_nDim * 2];
	for (int j = 0; j < 2 * m_ptrBEMEle->m_nDim; j++)
		m_dpRefCor[j] = 0;
	//generate reference node for torso surface
	for (int m = m_nNodNumInner; m < m_nNodNum; m++)
	{
		for (int n = 0; n < m_ptrBEMEle->m_nDim; n++)
			m_dpRefCor[n] += m_dpNodCor[3 * m + n];
	}
	for (int q = 0; q < m_ptrBEMEle->m_nDim; q++)
		m_dpRefCor[q] /= m_nNodNumOuter;

	//generate reference node for peri surface
	for (int m = 0; m < m_nNodNumInner; m++)
	{
		for (int n = 0; n < m_ptrBEMEle->m_nDim; n++)
			m_dpRefCor[n + m_ptrBEMEle->m_nDim] += m_dpNodCor[3 * m + n];
	}
	for (int q = 0; q < m_ptrBEMEle->m_nDim; q++)
		m_dpRefCor[q + m_ptrBEMEle->m_nDim] /= m_nNodNumInner;

}

BEM_nosource::~BEM_nosource()
{
}

//*******************************************************
//isOutward: test the firstly computed Nor is outward or inward
//if element on the torso surface, easy, just use a ref point within torso
//if element on heart suface, use norm to generate a point, then use tet to test whether that point is in 
//*******************************************************
int BEM_nosource::isOutward()
{
	//torso element
	if (m_ptrBEMEle->m_nEleSeq > m_nEleNumInner)
	{
		double dx = m_ptrBEMEle->m_dpNodCor[0] - m_dpRefCor[0];
		double dy = m_ptrBEMEle->m_dpNodCor[1] - m_dpRefCor[1];
		double dz = m_ptrBEMEle->m_dpNodCor[2] - m_dpRefCor[2];
		if (dx * m_ptrBEMEle->m_dpNor[0] + dy * m_ptrBEMEle->m_dpNor[1]
				+ dz * m_ptrBEMEle->m_dpNor[2] > 0)
			return 1;
		else
			return 0;
	}
	else
	{
		double dx = m_ptrBEMEle->m_dpNodCor[0] - m_dpRefCor[3];
		double dy = m_ptrBEMEle->m_dpNodCor[1] - m_dpRefCor[4];
		double dz = m_ptrBEMEle->m_dpNodCor[2] - m_dpRefCor[5];
		if (dx * m_ptrBEMEle->m_dpNor[0] + dy * m_ptrBEMEle->m_dpNor[1]
				+ dz * m_ptrBEMEle->m_dpNor[2] < 0)
			return 1;
		else
			return 0;
		/*		varCor.x = m_ptrBEMEle->m_dpNodCor[0] - m_ptrBEMEle->m_dpNor[0]/NORM_SCALE;
		 varCor.y = m_ptrBEMEle->m_dpNodCor[1] - m_ptrBEMEle->m_dpNor[1]/NORM_SCALE;
		 varCor.z = m_ptrBEMEle->m_dpNodCor[2] - m_ptrBEMEle->m_dpNor[2]/NORM_SCALE;
		 if(m_ptrMFree->isout(varCor))
		 return 1;
		 else
		 return 0;
		 */}

}
//***********************************************************
void BEM_nosource::eleNor(int neleSeq)
{
	m_ptrBEMEle->norG();
	if (!isOutward())
	{
		for (int j = 0; j < m_ptrBEMEle->m_nDim; j++)
			m_ptrBEMEle->m_dpNor[j] = -m_ptrBEMEle->m_dpNor[j];
	}
	int pos = (neleSeq - 1) * m_ptrBEMEle->m_nDim;
	for (int i = 0; i < m_ptrBEMEle->m_nDim; i++)
		m_dpEleDer[pos++] = m_ptrBEMEle->m_dpNor[i];
	m_dpEleAre[neleSeq - 1] = absolute(m_ptrBEMEle->m_nG) / 2.0;
}

//*********************************************************
//virtual functions
//*********************************************************
void BEM_nosource::intPK_first(int nnodSeq, int neleSeq)
{
	eleNor(neleSeq);
	//generate single K,P
	m_ptrBEMEle->sinKP(*this, isBel(nnodSeq, neleSeq), 1);
//	cout<<m_ptrBEMEle->m_dpsinK[0]<<"  "<<m_ptrBEMEle->m_dpsinK[1]<<"  "<<m_ptrBEMEle->m_dpsinK[2]<<endl;
	// assemble into integral K,P
	int pos = (neleSeq - 1) * m_ptrBEMEle->m_nNodNum;
	int seq = 0;
//	wriFil(".\\heart836torso700\\output\\eleK.bin",m_ptrBEMEle->m_dpsinK, sizeof(double), 3, "a+b");
	for (int i = 0; i < m_ptrBEMEle->m_nNodNum; i++)
	{
		seq = m_npNodSeq[pos++];
//       wriFil(".\\heart836torso700\\output\\ele.bin",&seq, sizeof(int), 1, "a+b");
		m_dpK[(nnodSeq - 1) * m_nNodNum + seq - 1] += m_ptrBEMEle->m_dpsinK[i];
		m_dpP[(nnodSeq - 1) * m_nNodNum + seq - 1] += m_ptrBEMEle->m_dpsinP[i];
	}

//	double* ptr = mxGetPr(m_dpK.GetData());
//	wriFil(".\\heart836torso700\\test\\testK.bin",ptr,sizeof(double),(m_nNodNum+1)*m_nNodNum,"wb");
}

//*******************************************************
//BEMCompute
//final interface with external application, generate final P,K,B,A and Transfer matrix
//*******************************************************
void BEM_nosource::BEMCompute(char** filOutput)
{
	clock_t s, f;
	double d;
	s = clock();
	//loop for each surface node
	cout << "********BEM*****************************************" << endl;
	for (int i = 1; i <= m_nNodNum; i++)
	{
		//loop for each element in surface
		if (!isForbidden(i))
		{
			cerr << i;
			m_ptrBEMEle->setNodSeq(i, *this);
			if (m_bisFirst)
			{
				for (int j = 1; j <= m_nEleNum; j++)
				{
					m_ptrBEMEle->setEleSeq(j, *this);
					intPK_first(i, j);
				}
				m_bisFirst = false;
			}
			else
			{
				for (int j = 1; j <= m_nEleNum; j++)
				{
					m_ptrBEMEle->setEleSeq(j, *this);
					intPK(i, j);
				}
			}
//		test("K_step.bin",&m_dpK[(i-1)*m_nNodNum],sizeof(double),m_nNodNum,"wb");
		}
	}
	//write out the "outward normal", to test whether it is really outward
	test("BEM_outward_normal.bin", m_dpEleDer, sizeof(double),
			m_nEleNum * m_ptrBEMEle->m_nDim, "wb");
	test("BEM_gp.cor", m_dpGps, sizeof(double), m_nGPNum * m_ptrBEMEle->m_nDim,
			"wb");
	test("original_BEMK.bin", m_dpK, sizeof(double), m_nNodNum * m_nNodNum,
			"wb");
	test("original_BEMP.bin", m_dpP, sizeof(double), m_nNodNum * m_nNodNum,
			"wb");

	rearrangePK();
	// get Hii from row sum elimination
	int k, pos = 0;
	int num = m_nNodNum - m_nForNum;
	double temp = 0;
	for (k = 0; k < num; k++)
	{
		temp = 0;
		pos = k * num;
		for (int l = 0; l < num; l++)
		{
			if (l != k)
				temp += m_dpK_arranged[pos + l];
		}
		m_dpK_arranged[pos + k] = -temp;
	}

//	mwIndex col;
//	m_dpKInner = m_dpK(col,colon(1,m_nNodNumInner);
//	m_dpKOuter = m_dpK(col,colon(m_nNodNumInner+1,m_nNodNum);
//	m_dpPInner = m_dpP(col,colon(1,m_nNodNumInner);

	//write out
	writeFile(filOutput[0], m_dpK_arranged, sizeof(double), num * num, "wb");
	writeFile(filOutput[1], m_dpP_arranged, sizeof(double), num * num, "wb");
	f = clock();
	d = double(f - s) / CLOCKS_PER_SEC;
	::cout << "time for BEM" << d << endl;
	//MFree compute B

}
