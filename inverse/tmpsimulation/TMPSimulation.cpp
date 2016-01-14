#include "inclu.h"
#include "external.h"
#include "matEg.h"
#include "Matrix.h"
#include "TMPSimulation.h"

extern matEg engine;

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

//****************************************************************************************************************
//class applicable for Poly (p,p_dx,p_dy,p_dz); Shape; A, B
//***************************************************************************************************************
//without dimension information, intialize as empty matrix
MatrixGroup::MatrixGroup() :
		d(), dx(), dy(), dz() {
}

MatrixGroup::MatrixGroup(int row, int col) :
		d(row, col), dx(row, col), dy(row, col), dz(row, col) {
}

MatrixGroup::~MatrixGroup() {
}

void MatrixGroup::Initialize(int row, int col) {
	d.setData(0.0, row, col);
	dx.setData(0.0, row, col);
	dy.setData(0.0, row, col);
	dz.setData(0.0, row, col);
}

void MatrixGroup::reSet() {
	d.reSet(0.0);
	dx.reSet(0.0);
	dy.reSet(0.0);
	dz.reSet(0.0);
}

void MatrixGroup::setData(double data, int row, int col) {
	d.setData(data, row, col);
	dx.setData(data, row, col);
	dy.setData(data, row, col);
	dz.setData(data, row, col);
}

//**************************************************************************************************************
//Definition of classes
//class: bgMesh
//**************************************************************************************************************
//*******************************************************
//ctors
//set the domain of the mesh as the first mesh , and set natural gp info for later use
// set the initial bgmesh as outside the original domain(i.e. isOut = true)
//*******************************************************
bgMesh::bgMesh(double* cor, int num) :
		m_dMesInt(BG_INTERVAL) {
	double min_x = cor[0], min_y = cor[1], min_z = cor[2], max_x = cor[0],
			max_y = cor[1], max_z = cor[2], dx = 0, dy = 0, dz = 0;
	int count = 0;
	for (int i = 0; i < num; i++) {
		if (cor[count] < min_x)
			min_x = cor[count];
		else if (cor[count] > max_x)
			max_x = cor[count];
		count++;
		if (cor[count] < min_y)
			min_y = cor[count];
		else if (cor[count] > max_y)
			max_y = cor[count];
		count++;
		if (cor[count] < min_z)
			min_z = cor[count];
		else if (cor[count] > max_z)
			max_z = cor[count];
		count++;
	}

	dx = max_x - min_x;
	dy = max_y - min_y;
	dz = max_z - min_z;

	m_nNumX = ceil(dx / m_dMesInt);
	m_nNumY = ceil(dy / m_dMesInt);
	m_nNumZ = ceil(dz / m_dMesInt);

	cout << "Number of mesh in x direction" << m_nNumX << "\n"
			<< "Number of mesh in y direction" << m_nNumY
			<< "\n" "Number of mesh in z direction" << m_nNumZ << endl;

	m_dTolX = (m_nNumX * m_dMesInt - dx) / 2.0;
	m_dTolY = (m_nNumY * m_dMesInt - dy) / 2.0;
	m_dTolZ = (m_nNumZ * m_dMesInt - dz) / 2.0;

	m_dgxl = min_x - m_dTolX;
	m_dgxu = max_x + m_dTolX;
	m_dgyl = min_y - m_dTolY;
	m_dgyu = max_y + m_dTolY;
	m_dgzl = min_z - m_dTolZ;
	m_dgzu = max_z + m_dTolZ;

	m_dxl = m_dxu = m_dgxl;
	m_dyl = m_dyu = m_dgyl;
	m_dzl = m_dzu = m_dgzl;

	m_dpNatGPCor = new double[MAX_NAT_GP];
	m_dpNatGPWei = new double[MAX_NAT_GP];
	m_dpGPCor = new double[MAX_GP * DIM];
	m_dpGPWei = new double[MAX_GP * DIM];

	m_nSta = 2;
	m_nGPNum = 0;
	m_nGPNumAll = 0;
	m_ngi = 0;
	m_ngj = 0;
	m_ngk = 0;
}

//*******************************************************
//dtor
//*******************************************************
bgMesh::~bgMesh() {

	delete[] m_dpNatGPCor;
	delete[] m_dpNatGPWei;
	delete[] m_dpGPCor;
	delete[] m_dpGPWei;
}

//***************************************************
//generate all the gpts inside the heart
//***************************************************
void bgMesh::generateGPAll(double* cor, int* tet, int num, int num_tet) {
	double* gpcor = new double[DIM];
	double* gpwei = new double[DIM];
	int count = 0;

	for (int i = 0; i < m_nNumX; i++)
		for (int j = 0; j < m_nNumY; j++)
			for (int k = 0; k < m_nNumZ; k++) {
				setMesSeq(i, j, k);
				comSta(cor, num);
				if (m_nSta != 0) //start to check gpts
						{
					for (int l = 0; l < m_nGPNum; l++)
						for (int m = 0; m < m_nGPNum; m++)
							for (int n = 0; n < m_nGPNum; n++) {
								genGloGP(gpcor, gpwei, l, m, n);
								if (!isOut(gpcor, cor, tet, num_tet)) {
									m_dpGPCor[count++] = gpcor[0];
									m_dpGPCor[count++] = gpcor[1];
									m_dpGPCor[count++] = gpcor[2];
									count -= 3;
									m_dpGPWei[count++] = gpwei[0];
									m_dpGPWei[count++] = gpwei[1];
									m_dpGPWei[count++] = gpwei[2];
									m_nGPNumAll++;
								}
							}
				}

			}

	delete[] gpcor;
	delete[] gpwei;

	if (count != 3 * m_nGPNumAll) {
		cout << "Error in generating gpts in bg_mesh!" << endl;
		exit(0);
	}

	::test("gp.cor", m_dpGPCor, sizeof(double), m_nGPNumAll * 3, "wb");
//	::test("gp.cor",m_dpGPCor,sizeof(double),m_nGPNumAll*3,"wb");
}

//*******************************************************
//setMesSeq
//change the mesh domain as loop continues, i.e. to sweep through all the meshes 
//*******************************************************
void bgMesh::setMesSeq(int i, int j, int k) {
	m_ngi = i;
	m_ngj = j;
	m_ngk = k;
	// generate domain of specified mesh
	m_dxl = m_dgxl + m_ngi * m_dMesInt;
	m_dxu = m_dxl + m_dMesInt;

	m_dyl = m_dgyl + m_ngj * m_dMesInt;
	m_dyu = m_dyl + m_dMesInt;

	m_dzl = m_dgzl + m_ngk * m_dMesInt;
	m_dzu = m_dzl + m_dMesInt;

	m_nSta = 2;
	m_nGPNum = 0;
}

//*******************************************************
//isout
//determine whether the mesh is totally inside the original domain
//set the "m_nSta", with the whole coordinate info from MFree class
// m_nSta = 0, while totally out
// =1, when partially out, in such case, gps' position needs to be checked
// else, totally in 
//*******************************************************
void bgMesh::comSta(double* cor, int num) {

	int count = 0;
	for (int i = 0; i < num; i++) {
		//at least one nodal is in the mesh, so it's not totally out
		if (cor[count] >= m_dgxl && cor[count] <= m_dgxu
				&& cor[count + 1] >= m_dgyl && cor[count + 1] <= m_dgyu
				&& cor[count + 2] >= m_dgzl && cor[count + 2] <= m_dgzu)
			m_nGPNum = m_nGPNum + 1;

		count += 3;
	}

	if (m_nGPNum == 0)
		m_nSta = 0;

	// generate natural gp coordinates&weights, just for 1D info

	else {
		//m_nGPNum = floor (pow(m_nGPNum, 1.0/3.0 ) + 2);
		m_nGPNum = NUM_GP;
		for (int j = 0; j < MAX_NAT_GP; j++) {
			m_dpNatGPCor[j] = 0;
			m_dpNatGPWei[j] = 0;
		}
		switch (m_nGPNum) {
		case 2:
			m_dpNatGPCor[0] = -0.577350269189626;
			m_dpNatGPCor[1] = 0.577350269189626;
			m_dpNatGPWei[0] = 1.0;
			m_dpNatGPWei[1] = 1.0;
			break;

		case 3:
			m_dpNatGPCor[0] = -0.774596669241483;
			m_dpNatGPCor[1] = 0;
			m_dpNatGPCor[2] = 0.774596669241483;
			m_dpNatGPWei[0] = 0.555555555555556;
			m_dpNatGPWei[1] = 0.888888888888889;
			m_dpNatGPWei[2] = 0.555555555555556;
			break;

		case 4:
			m_dpNatGPCor[0] = -0.861136311594053;
			m_dpNatGPCor[1] = -0.339981043584856;
			m_dpNatGPCor[2] = 0.339981043584856;
			m_dpNatGPCor[3] = 0.861136311594053;
			m_dpNatGPWei[0] = 0.347854845127454;
			m_dpNatGPWei[1] = 0.652145154862546;
			m_dpNatGPWei[2] = 0.347854845127454;
			m_dpNatGPWei[3] = 0.652145154862546;
			break;
		default:
			m_dpNatGPCor[0] = -0.906179845938664;
			m_dpNatGPCor[1] = -0.538469310105683;
			m_dpNatGPCor[2] = 0;
			m_dpNatGPCor[3] = 0.538469310105683;
			m_dpNatGPCor[4] = 0.906179845938664;
			m_dpNatGPWei[0] = 0.236926885056189;
			m_dpNatGPWei[1] = 0.478628670499366;
			m_dpNatGPWei[2] = 0.568888888888889;
			m_dpNatGPWei[3] = 0.478628670499366;
			m_dpNatGPWei[4] = 0.236926885056189;
			m_nGPNum = 5;
			break;
		}
	}

}

//*******************************************************
//genGloGP
//generate global gp coordinates and weights, one point each time,
//loop through all the gps in a mesh
//*******************************************************
void bgMesh::genGloGP(double* cor, double* wei, int l, int m, int n) {
	// x = ( xl + xu) /2 + gp*(xu-xl)/2
	cor[0] = (m_dxl + m_dxu) / 2.0 + m_dpNatGPCor[l] * (m_dxu - m_dxl) / 2.0;
	cor[1] = (m_dyl + m_dyu) / 2.0 + m_dpNatGPCor[m] * (m_dyu - m_dyl) / 2.0;
	cor[2] = (m_dzl + m_dzu) / 2.0 + m_dpNatGPCor[n] * (m_dzu - m_dzl) / 2.0;
	// neww = (xu-xl)*oldw/2
	wei[0] = m_dpNatGPWei[l] * (m_dxu - m_dxl) / 2.0;
	wei[1] = m_dpNatGPWei[m] * (m_dyu - m_dyl) / 2.0;
	wei[2] = m_dpNatGPWei[n] * (m_dzu - m_dzl) / 2.0;

}

//************************************************************
//isOut
//judge whether the gpts is in or outside of the tet
//************************************************************
bool bgMesh::isOut(double* gpcor, double* cor, int* tet, int num_tet) {
	//loop through all tetras, to find whether at least one include this gp
	int pos = 0;
	double* c0 = NULL;
	double* c1 = NULL;
	double* c2 = NULL;
	double* c3 = NULL;
	for (int i = 0; i < num_tet; i++) {
		int n0 = tet[pos++];
		int n1 = tet[pos++];
		int n2 = tet[pos++];
		int n3 = tet[pos++];

		c0 = &cor[(n0 - 1) * DIM];
		c1 = &cor[(n1 - 1) * DIM];
		c2 = &cor[(n2 - 1) * DIM];
		c3 = &cor[(n3 - 1) * DIM];

		int V0 = inVol(c0, c1, c2, c3);
		int V1 = inVol(gpcor, c1, c2, c3);
		int V2 = inVol(c0, gpcor, c2, c3);
		int V3 = inVol(c0, c1, gpcor, c3);
		int V4 = inVol(c0, c1, c2, gpcor);

		if (V0 == 0) {
			cerr << "The tetrahedron is degenerate !!" << endl;
			exit(0);
		} else {
			if (V0 > 0) {
				if (V1 >= 0 && V2 >= 0 && V3 >= 0 && V4 >= 0)
					return false; // true only if all are same sign as V0
			}	//else
				//	return true;
			else {
				if (V1 <= 0 && V2 <= 0 && V3 <= 0 && V4 <= 0)
					return false; // true only if all are same sign as V0
			}	//	else
			//		return true;
		}
	}
	return true;
}

//**************************************************************************************************************
//Definition of classes
//class: node
//**************************************************************************************************************
//*********************************************************
//ctors / dtor
//default ctor is needed for generating an array
//*********************************************************
node::node() {
	m_dpCor = new double[DIM];
	m_dpFibOrWei = new double[DIM];
	for (int i = 0; i < DIM; i++)
		m_dpCor[i] = 0;
	for (int i = 0; i < DIM; i++)
		m_dpFibOrWei[i] = 0;

	m_dD = DI_L;
}

node::node(double* cor, double* fib, double D) {
	m_dpCor = new double[DIM];
	m_dpFibOrWei = new double[DIM];
	for (int i = 0; i < DIM; i++)
		m_dpCor[i] = cor[i];
	for (int i = 0; i < DIM; i++)
		m_dpFibOrWei[i] = fib[i];

	m_dD = D;
}

node::~node() {
	delete[] m_dpCor;
	delete[] m_dpFibOrWei;
}

//**************************************************
// reInit: change the value of nodes
//***************************************************
void node::reInit(double* cor, double* fib, double D) {
	for (int i = 0; i < DIM; i++)
		m_dpCor[i] = cor[i];
	for (int i = 0; i < DIM; i++)
		m_dpFibOrWei[i] = fib[i];

	m_dD = D;
}

//*********************************************************************************************************************
//Definition of class
//class:nodeMF
//*********************************************************************************************************************
//**********************************
//ctors/dtor
//***********************************
nodeMF::nodeMF() :
		node(), m_mp(4, 1) {
	m_dRInf = 0;

	m_npIInf = new int[NUM_INF];
	for (int i = 0; i < NUM_INF; i++)
		m_npIInf[i] = 0;
}

nodeMF::nodeMF(double* cor, double* fib, double D, double* gcor, int num) :
		node(cor, fib, D), m_mp(4, 1) {
	m_npIInf = new int[NUM_INF];
	for (int i = 0; i < NUM_INF; i++)
		m_npIInf[i] = 0;
	calculRInf(gcor, num);
}

nodeMF::~nodeMF() {
	delete[] m_npIInf;
}

//***************************************
//reInit
//useful in initializng each node in an array
//****************************************
void nodeMF::reInit(double* cor, double* fib, double D, double* gcor, int num) {
	node::reInit(cor, fib, D);
	calculRInf(gcor, num);
	m_mp.setData(1.0, 0, 0);
	for (int i = 0; i < DIM; i++)
		m_mp.setData(m_dpCor[i], i + 1, 0);
}

//*************************************
//calculRInf
//Whenever a new node is specified, calculate its radius of inf domain
//*************************************
void nodeMF::calculRInf(double* cor, int num) {
	double r = 0;
	double* inf = new double[NUM_INF];
	for (int i = 0; i < NUM_INF; i++)
		inf[i] = MAX_R;
	int count = 0;

	//loop through every nodal
	for (int i = 0; i < num; i++) {
		double dx = cor[count++] - m_dpCor[0];
		double dy = cor[count++] - m_dpCor[1];
		double dz = cor[count++] - m_dpCor[2];
		r = sqrt(dx * dx + dy * dy + dz * dz);
		//check through the record queue, and add/update accordingly
		//determine the pos in record queue
		if (r != 0.0) {
			int k = 0;
			for (k = 0; k < NUM_INF && r > inf[k]; k++) {
			}
			//update record queue( insert in new neighbor at correct pos, and move following neighbors backward)
			if (k < NUM_INF) {
				for (int m = (NUM_INF - 1); m > k; m--) {
					inf[m] = inf[m - 1];
					m_npIInf[m] = m_npIInf[m - 1];
				}

				inf[k] = r;
				m_npIInf[k] = i;
			}
		}
	}

	m_dRInf = 0;
	for (int i = 0; i < NUM_INF; i++)
		m_dRInf += inf[i];
	m_dRInf = m_dRInf / NUM_INF;
//	m_dRInf = inf[NUM_INF-1];

	delete[] inf;
}

//*********************************************************************************************************************
//Definition of class
//class:nodeGS
//*********************************************************************************************************************
//**********************************
//ctors/dtor
//**********************************
//used in initializing an array of nodes, work with reInit()
nodeGS::nodeGS() :
		node(), m_mSha(), m_mDB(), m_mP(4, 1), m_mA(4, 4), m_mB(), m_mInvA(4,
				4), m_mD(3, 3), m_mQ(3, 3) {
	m_mQ.setData(1.0, 0, 0);
	m_mQ.setData(1.0, 1, 1);
	m_mQ.setData(1.0, 2, 2);

	if (IS_ISO == false) {
		m_mD.setData(DI_L, 0, 0);
		m_mD.setData(DI_S, 1, 1);
		m_mD.setData(DI_S, 2, 2);
	} else {
		m_mD.setData(DI, 0, 0);
		m_mD.setData(DI, 1, 1);
		m_mD.setData(DI, 2, 2);
	}

	m_dpFib = new double[DIM];
	for (int i = 0; i < DIM; i++)
		m_dpFib[i] = 0;
	m_nSur = 0;
	m_npISur = new int[MAX_INF];
	for (int i = 0; i < MAX_INF; i++)
		m_npISur[i] = 0;
	m_dpWSur = new double[MAX_INF * 4];
	for (int i = 0; i < 4 * MAX_INF; i++)
		m_dpWSur[i] = 0;
}

nodeGS::nodeGS(double* cor, double* wei, nodeMF* mfree, int num) :
		node(cor, wei, 0), m_mSha(), m_mDB(), m_mP(4, 1), m_mA(4, 4), m_mB(), m_mInvA(
				4, 4), m_mD(3, 3), m_mQ(3, 3) {
	m_mQ.setData(1.0, 0, 0);
	m_mQ.setData(1.0, 1, 1);
	m_mQ.setData(1.0, 2, 2);

	m_dpFib = new double[DIM];
	for (int i = 0; i < DIM; i++)
		m_dpFib[i] = 0;
	m_nSur = 0;
	m_npISur = new int[MAX_INF];
	for (int i = 0; i < MAX_INF; i++)
		m_npISur[i] = 0;
	m_dpWSur = new double[MAX_INF * 4];
	for (int i = 0; i < 4 * MAX_INF; i++)
		m_dpWSur[i] = 0;

	Poly3D(m_mP.d, cor[0], cor[1], cor[2]);
	Poly3D_dx(m_mP.dx, cor[0], cor[1], cor[2]);
	Poly3D_dy(m_mP.dy, cor[0], cor[1], cor[2]);
	Poly3D_dz(m_mP.dz, cor[0], cor[1], cor[2]);

	calInf(mfree, num);
	initMatrixInf();
	calSha(mfree);
	calFib(mfree);
	calRot();
	calDB(mfree);
}

nodeGS::~nodeGS() {
	delete[] m_dpFib;
	delete[] m_npISur;
	delete[] m_dpWSur;
}

//******************************************
//reInit:useful for initializing gpts in an array
//******************************************
void nodeGS::reInit(double* cor, double* wei, nodeMF* mfree, int num) {
	node::reInit(cor, wei, 0);
	Poly3D(m_mP.d, cor[0], cor[1], cor[2]);
	Poly3D_dx(m_mP.dx, cor[0], cor[1], cor[2]);
	Poly3D_dy(m_mP.dy, cor[0], cor[1], cor[2]);
	Poly3D_dz(m_mP.dz, cor[0], cor[1], cor[2]);

	calInf(mfree, num);
	initMatrixInf();
	calSha(mfree);
	calFib(mfree);
//	calRot();
	calRotbyRotation();
	calDB(mfree);
}

//*****************************************
//initMatrixInf: allocate the spaces for those matrices whose sizes are related to #supporting nodes
//*****************************************
void nodeGS::initMatrixInf() {
	m_mB.Initialize(4, m_nSur);
	m_mSha.setData(0.0, 4, m_nSur);
	m_mDB.setData(0.0, 3, m_nSur);
}

//*****************************************
//calInf: calculate support nodes & shape function for each gp
//******************************************
void nodeGS::calInf(nodeMF* mfree, int num) {
	double a, r, dx, dy, dz;
	int count = 0;
	//loop through each nodal, to check which ones influence domain cover gp
	for (int i = 0; i < num; i++) { //the distance between gp and current nodals
		a = DMAX * mfree[i].m_dRInf;
		dx = mfree[i].m_dpCor[0] - m_dpCor[0];
		dy = mfree[i].m_dpCor[1] - m_dpCor[1];
		dz = mfree[i].m_dpCor[2] - m_dpCor[2];
		r = sqrt(dx * dx + dy * dy + dz * dz);

		if (r > 1e-08) {
			dx /= (a * r);
			dy /= (a * r);
			dz /= (a * r);
			r /= a;

			if (r <= 0.5) {
				if (m_nSur + 1 > MAX_INF) {
					cout << "influencing nodes out of limits" << endl;
					exit(0);
				} else {
					// w = 2/3 - 4r^2 + 4r^3
					m_dpWSur[count++] = 2.0 / 3.0 - 4 * r * r + 4 * r * r * r;
					// wx = -8rx + 12r*rx    ( rx = (x - xi)/(a^2*r))
					m_dpWSur[count++] = (-8 * r + 12 * r * r) * dx;
					m_dpWSur[count++] = (-8 * r + 12 * r * r) * dy;
					m_dpWSur[count++] = (-8 * r + 12 * r * r) * dz;

					m_npISur[m_nSur++] = i + 1;
				}
			} else if (0.5 < r && r <= 1) {
				if (m_nSur + 1 > MAX_INF) {
					cout << "influencing nodes out of limits" << endl;
					wait();
					exit(0);
				} else {
					// w = 4/3 - 4r + 4r^2 - 4/3r^3
					m_dpWSur[count++] = 4.0 / 3.0 - 4 * r + 4 * r * r
							- 4.0 / 3.0 * r * r * r;
					// wx = -4rx/r + 8rx - 4rrx
					m_dpWSur[count++] = (-4 + 8 * r - 4 * r * r) * dx;
					m_dpWSur[count++] = (-4 + 8 * r - 4 * r * r) * dy;
					m_dpWSur[count++] = (-4 + 8 * r - 4 * r * r) * dz;

					m_npISur[m_nSur++] = i + 1;
				}
			}
		}
	}
	if (count != 4 * m_nSur) {
		cout << "Errors in recording supporting nodes weight for gpts!" << endl;
		exit(0);
	}

	if (m_nSur > MAX_INF) {
		cout << "Error: # supporting nodes excess MAX_INF" << endl;
		wait();
	}
}

//***********************************************
//calSha: calculate shape function for each node
//***********************************************
void nodeGS::calSha(nodeMF* node) {
	Matrix tmp(4, 4);
	int count = 0;

	//loop through every supportive point
	//construct A,B, and their derivatives
	for (int i = 0; i < m_nSur; i++) {
		//compute A,B and their derivatives nodal by nodal

		tmp = node[m_npISur[i] - 1].m_mp * (~(node[m_npISur[i] - 1].m_mp));
		m_mA.d += tmp * m_dpWSur[count++];
		m_mA.dx += tmp * m_dpWSur[count++];
		m_mA.dy += tmp * m_dpWSur[count++];
		m_mA.dz += tmp * m_dpWSur[count++];

		count -= 4;
		m_mB.d.setSub(node[m_npISur[i] - 1].m_mp * m_dpWSur[count++], 0, 3, i,
				i);
		m_mB.dx.setSub(node[m_npISur[i] - 1].m_mp * m_dpWSur[count++], 0, 3, i,
				i);
		m_mB.dy.setSub(node[m_npISur[i] - 1].m_mp * m_dpWSur[count++], 0, 3, i,
				i);
		m_mB.dz.setSub(node[m_npISur[i] - 1].m_mp * m_dpWSur[count++], 0, 3, i,
				i);
	}
	//compute ddSha

	m_mInvA.d = m_mA.d.inv_gsl();

	m_mSha.setSub((~(m_mP.d)) * m_mInvA.d * m_mB.d, 0, 0, 0, m_nSur - 1);
//	m_mSha.setSub((~(m_mP.d))*m_mInvA.d*m_mB.d.sub(0,3,0,m_nSur-1), 0, 0, 0, m_nSur-1);
//	sha.d.save("../BEMForward/GP/sha.bin","a+b");
	if (checkNan(m_mSha)) {
		cout << "Nan happens @ sha" << endl;
		wait();
	}

	int dsum = m_mSha.sub(0, 0, 0, m_nSur - 1).sum() - 1.0;
	if (dsum > 0.00001 || dsum < -0.00001) {
		cout << dsum << endl;
		cout << m_mSha.sub(0, 0, 0, m_nSur - 1).sum() << endl;
		cout << "Error: sum of shape is not equal to 1 on gpt" << endl;
		wait();
	}

	m_mInvA.dx = -m_mInvA.d * m_mA.dx * m_mInvA.d;
	m_mInvA.dy = -m_mInvA.d * m_mA.dy * m_mInvA.d;
	m_mInvA.dz = -m_mInvA.d * m_mA.dz * m_mInvA.d;

	/*	if (checkNan(m_mInvA.d))
	 {
	 cout<<"Nan happens @ invA"<<endl;
	 wait();
	 }
	 if (checkNan(A.dx))
	 {
	 cout<<"Nan happens @ A_dx"<<endl;
	 wait();
	 }
	 */
//	m_mSha.setSub((~m_mP.dx)*m_mInvA.d*m_mB.d.sub(0,3,0,m_nSur-1) + (~m_mP.d)*m_mInvA.dx*m_mB.d.sub(0,3,0,m_nSur-1) + (~m_mP.d)*m_mInvA.d*m_mB.dx.sub(0,3,0,m_nSur-1),1,1,0,m_nSur-1);
//	m_mSha.setSub((~m_mP.dy)*m_mInvA.d*m_mB.d.sub(0,3,0,m_nSur-1) + (~m_mP.d)*m_mInvA.dy*m_mB.d.sub(0,3,0,m_nSur-1) + (~m_mP.d)*m_mInvA.d*m_mB.dy.sub(0,3,0,m_nSur-1),2,2,0,m_nSur-1);
//	m_mSha.setSub((~m_mP.dz)*m_mInvA.d*m_mB.d.sub(0,3,0,m_nSur-1)+ (~m_mP.d)*m_mInvA.dz*m_mB.d.sub(0,3,0,m_nSur-1) + (~m_mP.d)*m_mInvA.d*m_mB.dz.sub(0,3,0,m_nSur-1),3,3,0,m_nSur-1);
	m_mSha.setSub(
			(~m_mP.dx) * m_mInvA.d * m_mB.d + (~m_mP.d) * m_mInvA.dx * m_mB.d
					+ (~m_mP.d) * m_mInvA.d * m_mB.dx, 1, 1, 0, m_nSur - 1);
	m_mSha.setSub(
			(~m_mP.dy) * m_mInvA.d * m_mB.d + (~m_mP.d) * m_mInvA.dy * m_mB.d
					+ (~m_mP.d) * m_mInvA.d * m_mB.dy, 2, 2, 0, m_nSur - 1);
	m_mSha.setSub(
			(~m_mP.dz) * m_mInvA.d * m_mB.d + (~m_mP.d) * m_mInvA.dz * m_mB.d
					+ (~m_mP.d) * m_mInvA.d * m_mB.dz, 3, 3, 0, m_nSur - 1);

	if (checkNan(m_mSha)) {
		cout << "Nan happens @ sha_d" << endl;
		wait();
	}
	/*	if (checkNan(m_mSha.dy))
	 {
	 cout<<"Nan happens @ sha_dy"<<endl;
	 wait();
	 }
	 if (checkNan(m_mSha.dz))
	 {
	 cout<<"Nan happens @ sha_dz"<<endl;
	 wait();
	 }*/
}

//****************************************
//calFib
//interpolate fib direction from other supporting nodes
//*****************************************
void nodeGS::calFib(nodeMF* mfree) {
	/*
	 Matrix fSur_x(m_nSur,1);
	 Matrix fSur_y(m_nSur,1);
	 Matrix fSur_z(m_nSur,1);
	 for(int i=0; i< m_nSur; i++)
	 fSur_x.setData(mfree[m_npISur[i]-1].m_dpFibOrWei[0],i,0);
	 for(int i=0; i< m_nSur; i++)
	 fSur_y.setData(mfree[m_npISur[i]-1].m_dpFibOrWei[1],i,0);
	 for(int i=0; i< m_nSur; i++)
	 fSur_z.setData(mfree[m_npISur[i]-1].m_dpFibOrWei[2],i,0);

	 m_dpFib[0]= (m_mSha.sub(0,0,0,m_nSur-1)*fSur_x)(0,0);
	 m_dpFib[1] = (m_mSha.sub(0,0,0,m_nSur-1)*fSur_y)(0,0);
	 m_dpFib[2] = (m_mSha.sub(0,0,0,m_nSur-1)*fSur_z)(0,0);
	 */
	m_dpFib[0] = 0;
	m_dpFib[1] = 0;
	m_dpFib[2] = 0;

	Matrix sha = m_mSha.sub(0, 0, 0, m_nSur - 1);
	double* ptr = sha.getDataPtr();
	for (int i = 0; i < m_nSur; i++) {
		int idx = m_npISur[i] - 1;
		m_dpFib[0] += ptr[i] * mfree[idx].m_dpFibOrWei[0];
		m_dpFib[1] += ptr[i] * mfree[idx].m_dpFibOrWei[1];
		m_dpFib[2] += ptr[i] * mfree[idx].m_dpFibOrWei[2];
	}
	double r = sqrt(
			m_dpFib[0] * m_dpFib[0] + m_dpFib[1] * m_dpFib[1]
					+ m_dpFib[2] * m_dpFib[2]);
	m_dpFib[0] = m_dpFib[0] / r;
	m_dpFib[1] = m_dpFib[1] / r;
	m_dpFib[2] = m_dpFib[2] / r;
}

void nodeGS::calRotbyRotation() {
	//fiber direction becomes local z-axis
	double r = sqrt(m_dpFib[1] * m_dpFib[1] + m_dpFib[0] * m_dpFib[0]);
	double theta_y = 0, theta_z = 0;
	//rotation around y-axis
	//The rotation angle about y axis
	if (r != 0) {
		theta_y = atan(m_dpFib[2] / r);
	} else {
		if (m_dpFib[2] > 0) {
			theta_y = PI / 2;
		} else if (m_dpFib[2] < 0) {
			theta_y = -PI / 2;
		}
	}

	//the rotation angle about z axis
	if (m_dpFib[0] > 0) {
		theta_z = atan(m_dpFib[1] / m_dpFib[0]);
	} else if (m_dpFib[0] < 0) {
		theta_z = PI + atan(m_dpFib[1] / m_dpFib[0]);
	} else {
		if (m_dpFib[1] > 0)
			theta_z = PI / 2;
		else if (m_dpFib[1] < 0)
			theta_z = -PI / 2;
		else
			theta_z = 0;
	}

	//rotation matrix
	Matrix Ty(3, 3);
//	double *TyPtr=Ty.GetPtr();
	Ty.setData(cos(theta_y), 0, 0);
	Ty.setData(sin(theta_y), 2, 0);
	Ty.setData(1.0, 1, 1);
	Ty.setData(-sin(theta_y), 0, 2);
	Ty.setData(cos(theta_y), 2, 2);

	/*    double TM_y[]={
	 cos(theta_y),    0, -sin(theta_y),
	 0,    1,     0,
	 sin(theta_y),    0, cos(theta_y),};
	 memcpy( TyPtr, TM_y, 9 * sizeof(double) );
	 */
	Matrix Tz(3, 3);
	//    double *TzPtr=Tz.GetPtr();
	Tz.setData(cos(theta_z), 0, 0);
	Tz.setData(sin(theta_z), 1, 0);
	Tz.setData(-sin(theta_z), 0, 1);
	Tz.setData(cos(theta_z), 1, 1);
	Tz.setData(1.0, 2, 2);

	m_mQ = ~(Tz * Ty);
}
/*	    double TM_z[]={ cosh, -sinh,  0,
 sinh,   cosh,  0,
 0,      0,  1, };
 memcpy( TzPtr, TM_z, 9 * sizeof(double) );

 MyMat matglo2loc(3,3);

 matglo2loc=Tz*Ty; //transformation matrix from global to local

 m_mQ=matglo2loc.transpose();//transformation matrix from local to global
 */
/* fiber direction to x-axis
 //The rotation angle about y axis
 double r = sqrt(mdpFib*fx+fy*fy);
 if(r!=0)
 {
 pointSet.pts[i].fib_ang_v = atan(fz/r);
 }
 else
 {
 if(fz>0)
 {
 pointSet.pts[i].fib_ang_v =PI/2;
 }
 else if(fz<0)
 {
 pointSet.pts[i].fib_ang_v =-PI/2;
 }
 else
 {
 cout<<" error in GetFNSForPoints"<<endl;
 }
 }


 //the rotation angle about z axis
 if(fx>0)
 {
 pointSet.pts[i].fib_ang_h=atan(fy/fx);
 }
 else if(fx<0)
 {
 pointSet.pts[i].fib_ang_h=PI+atan(fy/fx);
 }
 else
 {
 if(fy>0)
 pointSet.pts[i].fib_ang_h=PI/2;
 else if(fy<0)
 pointSet.pts[i].fib_ang_h=-PI/2;
 else
 pointSet.pts[i].fib_ang_h=0;
 }

 }

 free( use );
 if( !canUse )
 {
 cout<< "DeformableObject::GetFNSForPoints: not all gauss points can be used."<<endl;
 exit(0);
 }
 */

//****************************************
//calRot
//calculate the rotation matrix needed to get global D
//*****************************************
void nodeGS::calRot() {
	double e1, e2, e3;
	if (m_dpFib[2] != 0 && m_dpFib[0] != 0) {
		//local x
		e1 = m_dpFib[0];
		e2 = m_dpFib[1];
		e3 = -(e1 * e1 + e2 * e2) / (m_dpFib[2]);
		if (m_dpFib[2] < 0) {
			e1 = -e1;
			e2 = -e2;
			e3 = -e3;
		}
		m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
		m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
		m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
		//for local y axis
		e1 = -m_dpFib[1] * m_dpFib[1] / m_dpFib[0];
		e2 = m_dpFib[1];
		e3 = 0;
		if (m_dpFib[1] * m_dpFib[0] < 0) {
			e1 = -e1;
			e2 = -e2;
			e3 = -e3;
		}
		m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
		m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
		m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
		//for local z axis (longitudal fiber)
		e1 = m_dpFib[0];
		e2 = m_dpFib[1];
		e3 = m_dpFib[2];
		m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
		m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
		m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 2);
	} else {
		//fiber direction in x-y plane
		if (m_dpFib[2] == 0 && m_dpFib[0] != 0) {
			e1 = 0;
			e2 = 0;
			e3 = -1;
			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
			//for local y axis
			e1 = -m_dpFib[1] * m_dpFib[1] / m_dpFib[0];
			e2 = m_dpFib[1];
			e3 = 0;
			if (m_dpFib[0] * m_dpFib[1] < 0) {
				e1 = -e1;
				e2 = -e2;
				e3 = -e3;
			}
			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
			//for local z axis (longitudal fiber)
			e1 = m_dpFib[0];
			e2 = m_dpFib[1];
			e3 = m_dpFib[2];
			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 2);
		}
		//fiber direction in y-z plane
		else if (m_dpFib[0] == 0 && m_dpFib[2] != 0) {
			e1 = m_dpFib[0];
			e2 = m_dpFib[1];
			e3 = -(m_dpFib[0] * m_dpFib[0] + m_dpFib[1] * m_dpFib[1])
					/ (m_dpFib[2]);
			if (m_dpFib[1] * m_dpFib[2] < 0) {
				e1 = -e1;
				e2 = -e2;
				e3 = -e3;
			}
			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
			//for local y axis
			e1 = -1;
			e2 = 0;
			e3 = 0;

			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
			//for local z axis (longitudal fiber)
			e1 = m_dpFib[0];
			e2 = m_dpFib[1];
			e3 = m_dpFib[2];
			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 1);
		} else {
			if (m_dpFib[1] > 0) {
				e1 = 0;
				e2 = 0;
				e3 = -1;
			} else {
				e1 = 0;
				e2 = 0;
				e3 = 1;
			}
			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
			//for local y axis
			e1 = -1;
			e2 = 0;
			e3 = 0;

			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
			//for local z axis (longitudal fiber)
			e1 = m_dpFib[0];
			e2 = m_dpFib[1];
			e3 = m_dpFib[2];
			m_mQ.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
			m_mQ.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
			m_mQ.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 2);
		}
	}

//	m_mQ = ~m_mQ;
}

//****************************************
//calDB
//calculate DB = sqrt(D0)QB, where localK = (DB)'DB
//*****************************************
void nodeGS::calDB(nodeMF* mfree) {
	if (IS_HOMO)
		m_mDB = m_mD * m_mQ * (m_mSha.sub(1, 3, 0, m_nSur - 1));
	else {
		double d = 0;
		for (int i = 0; i < m_nSur; i++)
			d += mfree[m_npISur[i] - 1].m_dD * m_mSha(0, i);
		m_mD.setData(d, 2, 2);
		d /= Ransio;
		m_mD.setData(d, 0, 0);
		m_mD.setData(d, 1, 1);

		m_mDB = m_mD * m_mQ * m_mSha.sub(1, 3, 0, m_nSur - 1);
	}
}

//*****************************************
//getWei
//return the weight of wx*wy*wz
//******************************************
double nodeGS::getWei() {
	return m_dpFibOrWei[0] * m_dpFibOrWei[1] * m_dpFibOrWei[2];
}

//*********************************************************************************************************************
//Definition of class
//class:Heart. Top-level class
//*********************************************************************************************************************
//**********************************
//ctors/dtor
//**********************************
Heart::Heart(char* path_in) {
//	path = new char[MAX_FILENAME];
//	strcpy(path,input);

	path = path_in;
	m_dpTetCor = NULL;
	m_npTet = NULL;

	readGeo();

	readMF();

	genBG();

	cout << "Number of gtps" << bg->m_nGPNumAll << endl;

	genGS();

	int count = m_nMFNum * m_nMFNum;
	m_dpM = new double[count];
	m_dpK = new double[count];
}

Heart::~Heart() {
//	delete []path;
	delete[] m_dpTetCor;
	delete[] m_npTet;
	delete[] m_dpM;
	delete[] m_dpK;

	delete[] mfree;
	delete[] gspt;
	delete bg;

}

//********************************
//readGeo
//read in information about tetra
//**********************************
void Heart::readGeo() {
	char* file = new char[MAX_FILENAME];
	strcpy(file, path);
	strcat(file, "Ref/heart_tet.cor");
	::readFile(file, &m_dpTetCor, sizeof(double), DIM, m_nTetVerNum);

	strcpy(file, path);
	strcat(file, "Ref/heart.tet");
	::readFile(file, &m_npTet, sizeof(int), 4, m_nTetNum);

	delete[] file;
}

//****************************************
//readMF
//read in mfree nodes, form array of mfree class
//****************************************
void Heart::readMF() {
	double* cor = NULL;
	double* fib = NULL;
	double* D = NULL;

	char sid[5];
	sprintf(sid, "%d", ID);
	char* file = new char[MAX_FILENAME];
	strcpy(file, path);
	strcat(file, "Final/");
	strcat(file, sid);
	strcat(file, "/heart.cor");
	::readFile(file, &cor, sizeof(double), DIM, m_nMFNum);
	/* 
    strcpy(file, path);
	strcat(file,"Final/");
	strcat(file,sid);
	strcat(file, "/heart.fib");
	::readFile(file,&fib, sizeof(double),DIM,m_nMFNum);
	*/
    if (!IS_HOMO) {
		strcpy(file, path);
		strcat(file, "Final/");
		strcat(file, sid);
		strcat(file, "/PJN.map");
		::readFile(file, &D, sizeof(double), 1, m_nMFNum);

		for (int i = 0; i < m_nMFNum; i++) {
			if (D[i] == 0)
				D[i] = DI_L;
			else
				D[i] = DI_L * P2M;
		}
	} else {
		D = new double[m_nMFNum];
		for (int i = 0; i < m_nMFNum; i++)
			D[i] = DI_L;
	}

	mfree = new nodeMF[m_nMFNum];
	double* ptrc = cor;
	double* ptrf = fib;

	for (int i = 0; i < m_nMFNum; i++) {
		mfree[i].reInit(ptrc, ptrc, D[i],cor, m_nMFNum);
		ptrc += 3;
		//ptrf+=3;
	}

	delete[] cor;
	delete[] fib;
	delete[] D;
	delete[] file;
}

//*******************************************
//genBG
//generate bg, and get all gs
//*******************************************
void Heart::genBG() {
	bg = new bgMesh(m_dpTetCor, m_nTetVerNum);
	bg->generateGPAll(m_dpTetCor, m_npTet, m_nTetVerNum, m_nTetNum);
}

//********************************************
//genGS
//generate gs, initialize
//********************************************
void Heart::genGS() {
	gspt = new nodeGS[bg->m_nGPNumAll];

	double* ptrc = bg->m_dpGPCor;
	double* ptrw = bg->m_dpGPWei;

	for (int i = 0; i < bg->m_nGPNumAll; i++) {
		gspt[i].reInit(ptrc, ptrw, mfree, m_nMFNum);
		ptrc += 3;
		ptrw += 3;
	}
}

//********************************************
//assTrans
//everything is ready, then assememble M, K
//********************************************
void Heart::assTrans() {
	// assemble M = Sha'Sha
	for (int i = 0; i < m_nMFNum * m_nMFNum; i++)
		m_dpM[i] = 0;
	for (int i = 0; i < bg->m_nGPNumAll; i++) {
		Matrix sha = gspt[i].m_mSha.sub(0, 0, 0, gspt[i].m_nSur - 1);
		double* ptr = sha.getDataPtr();	// gspt[i].m_mSha.getDataPtr();
		double wei = gspt[i].getWei();
		for (int j = 0; j < gspt[i].m_nSur; j++)
			for (int k = 0; k < gspt[i].m_nSur; k++)
				m_dpM[m_nMFNum * (gspt[i].m_npISur[k] - 1) + gspt[i].m_npISur[j]
						- 1] += wei * ptr[j] * ptr[k];	//wei*ptr[4*j]*ptr[4*k]
	}

	//assemble K
	for (int i = 0; i < m_nMFNum * m_nMFNum; i++)
		m_dpK[i] = 0;
	for (int i = 0; i < bg->m_nGPNumAll; i++) {
		double wei = gspt[i].getWei();
		Matrix tmp = ((~gspt[i].m_mDB) * gspt[i].m_mDB) * wei;
		double* ptr = tmp.getDataPtr();
		int count = 0;

		for (int j = 0; j < gspt[i].m_nSur; j++)
			for (int k = 0; k < gspt[i].m_nSur; k++)
				m_dpK[m_nMFNum * (gspt[i].m_npISur[j] - 1) + gspt[i].m_npISur[k]
						- 1] += ptr[count++];
	}
}

void Heart::writeTrans() {
	int count = m_nMFNum * m_nMFNum;

	cout << "writing files" << endl;
	char sid[5];
	sprintf(sid, "%d", ID);
	char* file = new char[MAX_FILENAME];
	strcpy(file, path);
	strcat(file, "Simulation/");
	strcat(file, sid);
	strcat(file, "/Input/Mass.bin");
	::writeFile(file, m_dpM, sizeof(double), count, "wb");

	strcpy(file, path);
	strcat(file, "Simulation/");
	strcat(file, sid);
	strcat(file, "/Input/Stiff.bin");
	::writeFile(file, m_dpK, sizeof(double), count, "wb");

	delete[] file;
}

//***********************************************
//output infR for each mfree node, for checking
//************************************************
void Heart::checkInfR(char* testPath) {
	double* data = new double[m_nMFNum];
	for (int i = 0; i < m_nMFNum; i++)
		data[i] = mfree[i].m_dRInf;

	::writeFile(testPath, data, sizeof(double), m_nMFNum, "wb");
}

//**************************************************
//output gp cor for checking
//***************************************************
void Heart::checkGP(char* testPath) {
	char* file = new char[MAX_FILENAME];
	//cor
	strcpy(file, testPath);
	strcat(file, "gpts.cor");
	::writeFile(file, bg->m_dpGPCor, sizeof(double), DIM * bg->m_nGPNumAll,
			"wb");

	//fib
	double* data = new double[bg->m_nGPNumAll * DIM];
	for (int i = 0; i < bg->m_nGPNumAll; i++) {
		data[DIM * i] = gspt[i].m_dpFib[0];
		data[DIM * i + 1] = gspt[i].m_dpFib[1];
		data[DIM * i + 2] = gspt[i].m_dpFib[2];
	}
	strcpy(file, testPath);
	strcat(file, "gpts.fib");
	::writeFile(file, data, sizeof(double), DIM * bg->m_nGPNumAll, "wb");
	//# of supporting nodes
	int* data1 = new int[bg->m_nGPNumAll];
	for (int i = 0; i < bg->m_nGPNumAll; i++)
		data1[i] = gspt[i].m_nSur;
	strcpy(file, testPath);
	strcat(file, "gpts.inf");
	::writeFile(file, data1, sizeof(int), bg->m_nGPNumAll, "wb");
	//seq of supporting nodes
	strcpy(file, testPath);
	strcat(file, "gpts_inf.seq");
	for (int i = 0; i < bg->m_nGPNumAll; i++)
		::writeFile(file, gspt[i].m_npISur, sizeof(int), gspt[i].m_nSur, "a+b");

	delete[] data;
	delete[] data1;
	delete[] file;
}

//*************************************
//getInf
//write the inf nodes idx, for testing
//*************************************
void Heart::checkInfMfree(char* path) {
	char* file = new char[MAX_FILENAME];

	strcpy(file, path);
	strcat(file, "mfree_inf.seq");

	for (int i = 0; i < m_nMFNum; i++)
		::writeFile(file, mfree[i].m_npIInf, sizeof(int), NUM_INF, "a+b");

	delete[] file;
}

//*********************************************************************************************************************
//Definition of class
//class:Simulator
//*********************************************************************************************************************
//**********************************
//ctors/dtor
//**********************************
Simulator::Simulator(char* path) :
		U(ID, 1), V(ID, 1), T(ID, ID), H(YDIM, ID), Sti_all(), Sti(ID, 1), Par(
				ID, 1), state(), mea(), u0_RK(ID, 1), v0_RK(ID, 1), ui_RK(ID,
				1), vi_RK(ID, 1), u0_fg(ID, 1), v0_fg(ID, 1), vec_one(ID, 1) {
	char sid[5];
	sprintf(sid, "%d", ID);
	char* in = new char[MAX_FILENAME];
	strcpy(in, path);
	strcat(in, "Simulation/");
	strcat(in, sid);
	strcat(in, "/Input/");
	//***read T***********************
	char* file = new char[MAX_FILENAME];

	char subfolder[MAX_FILENAME];
	printf("Enter the subfolder for Simulation input path: \n");
	scanf("%s", &subfolder);
	strcat(in, subfolder);

	cout << "Read Trans_state" << endl;
	double* data = NULL;
	int tmp;
	strcpy(file, in);
	strcat(file, "Trans_state.bin");
	::readFile(file, &data, sizeof(double), ID, tmp);
	if (tmp != ID) {
		cout << "Errors in reading Trans_state!" << endl;
		wait();
	}
	double* ptr = T.getDataPtr();
	for (int i = 0; i < ID * ID; i++)
		ptr[i] = data[i];
	delete[] data;
	data = NULL;
	engine.addData(T, "C");
	//***read H***********************
	cout << "Read Trans" << endl;
	strcpy(file, in);
	strcat(file, "Trans.bin");
	::readFile(file, &data, sizeof(double), ID, tmp);
	if (tmp != YDIM) {
		cout << "Errors in reading Trans!" << endl;
		wait();
	}
	ptr = H.getDataPtr();
	for (int i = 0; i < ID * YDIM; i++)
		ptr[i] = data[i];
	delete[] data;
	data = NULL;
	engine.addData(H, "H");
	//***read t***********************
	cout << "Read t? (y/n)" << endl;
	char c;
	cin >> c;
	if (c == 'y' || c == 'Y') {
		t = NULL;
		strcpy(file, in);
		strcat(file, "time");
		if (FHN_TYPE == 'P')
			strcat(file, "_P.bin");
		else if (FHN_TYPE == 'M')
			strcat(file, "_M.bin");
		else
			strcat(file, "_O.bin");
		readFile(file, &t, sizeof(double), 1, step_simul);
		step_save = step_simul;
		tmax = t[step_simul - 1];
	} else {
		t = NULL;
		cout << "Input the total time length" << endl;
		cin >> tmax;
		cout << "Input the temproal resolution for simulation" << endl;
		cin >> dt_simul;
		cout << "Input the temproal resolution for saving results" << endl;
		cin >> dt_save;
		step_simul = round(tmax / dt_simul) + 1;
		step_save = round(tmax / dt_save) + 1;
	}

	//***Initialize state & mea*************
	state.setData(0.0, ID, step_save);
	mea.setData(0.0, YDIM, step_save);
	Sti_all.setData(0.0, ID, step_save - 1);
	vec_one.reSet(1.0);
	sti_l = new int[step_simul - 1];
	sti_r = new int[step_simul - 1];
	sti_p = new int[step_simul - 1];
	for (int i = 0; i < step_simul - 1; i++) {
		sti_l[i] = 0;
		sti_r[i] = 0;
		sti_p[i] = 0;
	}

	cout << "Read excitation sites" << endl;
	setExcitation(path);

	cout << "Read pacing sites" << endl;
	setPacing(path);

	parameter = NULL;
	strcpy(file, in);
	cout << "Input parameter file name? (y/n)" << endl;
	cin >> c;
	if (c == 'y' || c == 'Y') {
		char filename[50];
		scanf("%s", filename);
		strcat(file, filename);
	} else
		strcat(file, "Parameter.bin");

	readFile(file, &parameter, sizeof(double), 1, tmp);
	if (tmp != ID) {
		cout << "Error in reading parameter!" << endl;
		wait();
	}
	Par.setData(parameter);

	delete[] file;
	delete[] in;

	//*********AT,RT*************************
	AT = new double[ID];
	for (int i = 0; i < ID; i++)
		AT[i] = 0;
	AT_idx = new int[ID];
	for (int i = 0; i < ID; i++)
		AT_idx[i] = 0;
	RT = new double[ID];
	for (int i = 0; i < ID; i++)
		RT[i] = 0;
	RT_idx = new int[ID];
	for (int i = 0; i < ID; i++)
		RT_idx[i] = 0;
	APD = new double[ID];
	for (int i = 0; i < ID; i++)
		APD[i] = 0;
}

Simulator::~Simulator() {
	delete[] t;
	delete[] exc;
	delete parameter;
	if (focus)
		delete[] focus;

	delete[] AT;
	delete[] AT_idx;
	delete[] RT;
	delete[] RT_idx;
	delete[] APD;

	delete[] sti_l;
	delete[] sti_r;
	delete[] sti_p;
}

//**********************************
//set normal excitation sites
//**********************************
void Simulator::setExcitation(char* path) {
	char sid[5];
	sprintf(sid, "%d", ID);
	//***read exc & generate Init_X according the scheme chosen***********************
	exc = new int[ID];
	num_exc = 0;
	char* file = new char[MAX_FILENAME];

	int* s = NULL;

	//read in excitation locations
	cout << "Normal excitation? (y/n)" << endl;
	char c;
	cin >> c;
	if (c == 'y' || c == 'Y') {
		strcpy(file, path);
		strcat(file, "Final/");
		strcat(file, sid);
		char subfolder[MAX_FILENAME];
		printf("Enter the file name for exc: \n");
		scanf("%s", &subfolder);
		strcat(file, "/");
		strcat(file, subfolder);
		int tmp = 0;
		readFile(file, &s, sizeof(int), 1, tmp);
		if (tmp != ID + 1) {
			cout << "Errors in reading exc!" << endl;
			wait();
		}
		num_exc = s[0];
		for (int i = 1; i < tmp; i++)
			exc[i - 1] = s[i];
		delete[] s;
	} else {
		for (int i = 0; i < ID; i++)
			exc[i] = 0;
	}

	//determine the starting & duration of stimulus
	step_delay = 0;
	cout << "input stimulus duration (ms)" << endl;
	cin >> duration_exc;
	cout
			<< "Input the delay of RV excitation relative to LV: (-1: RBBB; -2: LBBB; >=0, no block"
			<< endl;
	cin >> t_delay;
	if (t_delay == -2) {
		step_delay = -2;
//		genInitRV();
		processA('R');
		int i = 0;
		if (t) {
			while (duration_exc >= t[i]) {
				sti_r[i] = 1;
				i++;
			}

		} else {
			for (int i = 0; i <= round(duration_exc / (dt_simul)); i++)
				sti_r[i] = 1;
		}
	} else if (t_delay == -1) {
		step_delay = -1;
		//   	genInitLV();
		processA('L');
		int i = 0;
		if (t) {
			while (duration_exc >= t[i]) {
				sti_l[i] = 1;
				i++;
			}
		} else {
			for (int i = 0; i <= round(duration_exc / (dt_simul)); i++)
				sti_l[i] = 1;
		}
	} else {
		step_delay = 0;
		if (t)
			for (step_delay = 0;
					step_delay < step_simul && t[step_delay] < t_delay;
					step_delay++) {
			}
		else
			step_delay = round(t_delay / dt_simul);

		cout << "RV is delayed to LV by: " << step_delay << " steps \n" << endl;
		//  		genInitLV();
		processA('L');
		if (t) {
			int i = 0;
			while (duration_exc >= t[i]) {
				sti_l[i] = 1;
				i++;
			}
			i = step_delay;
			while (duration_exc + t_delay >= t[i]) {
				sti_r[i] = 1;
				i++;
			}
		} else {
			for (int i = 0; i <= round(duration_exc / dt_simul); i++)
				sti_l[i] = 1;
			for (int i = step_delay;
					i <= step_delay + round(duration_exc / dt_simul); i++)
				sti_r[i] = 1;
		}
	}
//	engine.addData(T,"C");	

	delete[] file;
}

//****************************************
//set pacing sites
//****************************************
void Simulator::setPacing(char* path) {
	char sid[5];
	sprintf(sid, "%d", ID);
	//***read exc & generate Init_X according the scheme chosen***********************
	char* file = new char[MAX_FILENAME];

	//***focus,par ***********************
	cout << "Input the beginning instant of pacing sites:" << endl;
	start_focus = 0;
	cin >> start_focus;
	for (int i = 0; i < step_simul; i++) {
		if (t) {
			if (start_focus >= t[i]) {
				startstep_focus = i;
				break;
			}
		} else {
			startstep_focus = round(start_focus / dt_simul);
		}
	}
	//read the pacing sites	
	focus = NULL;
	num_focus = 0;
	char c;
	cout << "Input pacing sites? (y/n)" << endl;
	cin >> c;
	if (c == 'y' || c == 'Y') {
		strcpy(file, path);
		strcat(file, "Final/");
		strcat(file, sid);
		cout << "input pacing sites file name" << endl;
		char name[20];
		scanf("%s", name);
		strcat(file, "/");
		strcat(file, name);
		int tmp = 0;
		readFile(file, &focus, sizeof(double), 1, tmp);
		if (tmp != ID) {
			cout << "Reading pacing sites error!" << endl;
			wait();
		}
		for (int i = 0; i < ID; i++) {
			if (focus[i] == 1)
				num_focus++;
		}
		cout << "number of pacing sites: " << num_focus << endl;

//		genPacing();
//		processA('P');

		cout << "Input pacing frequency" << endl;
		cin >> freq_focus;
		cout << "Input pacing duration" << endl;
		cin >> duration_focus;

		int pos = startstep_focus;

		if (t) {
			while (t[pos] <= start_focus + duration_focus) {
				sti_p[pos] = 1;
				pos++;
			}
		} else {
			while (pos < step_simul - 1) {
				if (pos + round(duration_focus / dt_simul) >= step_simul - 2) {
					for (int i = pos; i < step_simul - 1; i++)
						sti_p[i] = 1;
				} else {
					for (int i = pos;
							i <= pos + round(duration_focus / dt_simul); i++)
						sti_p[i] = 1;
				}
				pos = pos + round(duration_focus / dt_simul)
						+ round((1.0 / freq_focus * 1000) / dt_simul);
			}
		}
	}

//	engine.addData(T,"C");	

	delete[] file;
}

//***************************************************
//given seq of source nodes, generate Init_state in vector form
//***************************************************
void Simulator::genInitLV() {
	cout << "# of pacemakers: " << num_exc << endl;
	num_exc_lv = 0;
	if (SCHEME == 2) {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 1 || exc[i] == 3) {
				Sti.setData(U_INIT, i, 0);
				num_exc_lv++;
			}
		}
	} else {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 1 || exc[i] == 3) {
				state.setData(U_INIT, i, 0);
				num_exc_lv++;
			}
		}
	}
	cout << "# of pacemakers in LV: " << num_exc_lv << endl;
}

//***************************************************
//given seq of source nodes, for the delayed rv part
//***************************************************
void Simulator::genInitRV() {
	cout << "# of pacemakers: " << num_exc << endl;
	num_exc_rv = 0;
	if (SCHEME == 2) {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 2 || exc[i] == 4) {
				Sti.setData(U_INIT, i, 0);
				num_exc_rv++;
			}
		}
	} else {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 2 || exc[i] == 4) {
				state.setData(U_INIT, i, 0);
				num_exc_rv++;
			}
		}
	}
	cout << "# of pacemakers in RV: " << num_exc_rv << endl;
}

//***************************************************
//given seq of pacing sites, generate the pacing Stimulus vector
//***************************************************
void Simulator::genPacing() {
	Sti.reSet(0.0);

	for (int i = 0; i < ID; i++) {
		if (focus[i] == 1.0)
			Sti.setData(U_INIT, i, 0);
	}
}

//***************************************************
//processA
//assign rows corresponding to source node as zero
//***************************************************
void Simulator::processA(char lr) {
	if (lr == 'L') {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 1 || exc[i] == 3)
				T.setSub(Matrix(1, ID), i, i, 0, ID - 1);
		}
	} else if (lr == 'R') {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 2 || exc[i] == 4)
				T.setSub(Matrix(1, ID), i, i, 0, ID - 1);
		}
	} else {
		for (int i = 0; i < ID; i++) {
			if (focus[i] == 1)
				T.setSub(Matrix(1, ID), i, i, 0, ID - 1);
		}
	}

	engine.addData(T, "C");
}

//*************************************************
//removeSti
//remove the stimuli after certain time length
//*************************************************
void Simulator::removeSti(char lr) {
	if (lr == 'L') {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 1 || exc[i] == 3)
				Sti.setData(0.0, i, 0);
		}
	}
	if (lr == 'R') {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 2 || exc[i] == 4)
				Sti.setData(0.0, i, 0);
		}
	}
	if (lr == 'F') {
		for (int i = 0; i < num_focus; i++) {
			Sti.setData(0.0, focus[i], 0);
		}
	}
}

//*************************************************
//addSti
//add the stimuli after certain time length
//*************************************************
void Simulator::addSti(char lr) {
	if (lr == 'L') {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 1 || exc[i] == 3)
				Sti.setData(U_INIT, i, 0);
		}
	}
	if (lr == 'R') {
		for (int i = 0; i < ID; i++) {
			if (exc[i] == 2 || exc[i] == 4)
				Sti.setData(U_INIT, i, 0);
		}
	}
	if (lr == 'F') {
		for (int i = 0; i < num_focus; i++) {
			Sti.setData(U_INIT, focus[i], 0);
		}
	}
}

//***************************************************
// from sti_l/_r/_p, set the stimulus vector at time instant ins
//***************************************************
void Simulator::getCurrentSti(int ins) {
	Sti.reSet(0.0);

	if (sti_l[ins] == 1) {
		for (int i = 0; i < ID; i++)
			if (exc[i] == 1 || exc[i] == 3)
				Sti.setData(U_INIT, i, 0);
	}
	if (sti_r[ins] == 1) {
		for (int i = 0; i < ID; i++)
			if (exc[i] == 2 || exc[i] == 4)
				Sti.setData(U_INIT, i, 0);
	}
	if (sti_p[ins] == 1) {
		for (int i = 0; i < ID; i++)
			if (focus[i] == 1)
				Sti.setData(U_INIT, i, 0);
	}
}

//***************************************************
//
//**************************************************
void Simulator::updateA(int ins) {
	if (step_delay >= 0 && ins == step_delay)
		processA('R');

	if (num_focus > 0 && ins == startstep_focus)
		processA('P');
}

//***************************************************
//f & g
// ode for FHN (non-uniform prameter)
//***************************************************
void Simulator::fg(Matrix& u, Matrix& v) {
	u0_fg = u;
	v0_fg = v;

	if (PAR_TYPE == 'a') {
		if (FHN_TYPE == 'O') {
			u = (u0_fg & (vec_one - u0_fg) & (u0_fg - Par)) - v0_fg + Sti;
			v = (u0_fg - v0_fg * PAR_DO) * PAR_BO;
		} else if (FHN_TYPE == 'M') {
			u = (u0_fg & (vec_one - u0_fg) & (u0_fg - Par)) * PAR_C1
					- (u0_fg & v0_fg) * PAR_C2 + Sti;
			v = (u0_fg - v0_fg * PAR_DM) * PAR_BM;
		} else {
			u = (u0_fg & (vec_one - u0_fg) & (u0_fg - Par)) * PAR_K
					- (u0_fg & v0_fg) + Sti;
			v = ((u0_fg & (u0_fg - Par - vec_one)) * PAR_K + v0_fg) * (-PAR_E);
		}

	} else {
		if (FHN_TYPE == 'O') {
			u = (u0_fg & (vec_one - u0_fg) & (u0_fg - vec_one * PAR_AO)) - v0_fg
					+ Sti;
			v = Par & (u0_fg - v0_fg * PAR_DO);
		} else if (FHN_TYPE == 'M') {
			u = (u0_fg & (vec_one - u0_fg) & (u0_fg - vec_one * PAR_AM))
					* PAR_C1 - (u0_fg & v0_fg) * PAR_C2 + Sti;
			v = Par & (u0_fg - v0_fg * PAR_DM);
		} else {
			u = (u0_fg & (vec_one - u0_fg) & (u0_fg - vec_one * Par)) * PAR_K
					- (u0_fg & v0_fg) + Sti;
			v = Par
					& ((u0_fg & (u0_fg - vec_one * (Par + vec_one))) * PAR_K
							+ v0_fg) * (-1);
		}
	}
	engine.addData(u0_fg, "u0");
	engine.evalString("u0 = C*u0;");
	engine.getData(u0_fg, "u0");
	u += u0_fg;
}

//*******************************************
//RK
// batch process
// first check the [i,j]=size(x); i determine whether joint, j determine whether batch (j=1, single)
//*******************************************
void Simulator::samPro_RK(int idx, double dt) {
	u0_RK = U;
	v0_RK = V;
	ui_RK = U;
	vi_RK = V;
	double a = 0;

	fg(u0_RK, v0_RK);
//	char sid[5];
//	sprintf(sid,"%d",idx);
//	u0_RK.checkNanInf(sid);

	a = dt / 6.0;
	U += u0_RK * a;
	V += v0_RK * a;
	a = dt / 2.0;
	u0_RK *= a;
	u0_RK += ui_RK;
	v0_RK *= a;
	v0_RK += vi_RK;
	fg(u0_RK, v0_RK);
//	u0_RK.checkNanInf(sid);

	a = dt / 3.0;
	U += u0_RK * a;
	V += v0_RK * a;
	a = dt / 2.0;
	u0_RK *= a;
	u0_RK += ui_RK;
	v0_RK *= a;
	v0_RK += vi_RK;
	fg(u0_RK, v0_RK);
//	u0_RK.checkNanInf(sid);

	a = dt / 3.0;
	U += u0_RK * a;
	V += v0_RK * a;
	u0_RK *= dt;
	u0_RK += ui_RK;
	v0_RK *= dt;
	v0_RK += vi_RK;
	fg(u0_RK, v0_RK);
//	u0_RK.checkNanInf(sid);

	a = dt / 6.0;
	U += u0_RK * a;
	V += v0_RK * a;

//	state.setSub(U,0,ID-1,idx,idx);
//	Sti_all.setSub(Sti,0,ID-1,idx-1,idx-1);

	if (idx % ((int) (dt_save / dt_simul)) == 0) {
//		cout<<"idx = "<<idx<<endl;
		int iidx = idx / ((int) (dt_save / dt_simul));
		state.setSub(U, 0, ID - 1, iidx, iidx);
		Sti_all.setSub(Sti, 0, ID - 1, iidx - 1, iidx - 1);
//		U.save("/Users/maomaowlw/Research/Mywork/Data/Geometry/MICCAI10/Processed/Case2/Simulation/1045/Output/repi/TMP.bin", "a+b");
//		state.setSub(v,ID,2*ID-1,iidx,iidx);
	}

}

//***************************************************
//staPro
//forward state model simulation with RK
//***************************************************
void Simulator::staPro() {
	for (int i = 1; i < step_simul; i++) {
		if (t) {
			dt_simul = t[i] - t[i - 1];
			dt_save = dt_simul;
		}

//		cout<<"Get Sti"<<endl;
		getCurrentSti(i - 1);
//		cout<<"update A"<<endl;
		updateA(i - 1);
//		cout<<"RK"<<endl;
		samPro_RK(i, dt_simul);
	}
	//else
	/*	else
	 {
	 for(int i=1; i<step; i++)
	 {
	 dt= t[i] - t[i-1];
	 if(num_focus>0 && i == startstep_focus+1)
	 processA('P');
	 if(i==step_delay+1)
	 processA('R');
	 samPro_RK(i,dt);
	 }
	 }*/
	//  	engine.addData(H,"H");
	engine.addData(state, "staX");
//   	engine.evalString("[r,c] = size(staX);");
	engine.evalString("Y=H*staX;");
	engine.getData(mea, "Y");
//	mea = H*staX.sub(0,UDIM-1,0,STEP-1);
}

//*********************************************************
//numDiff
//********************************************************
double Simulator::numDiff(double* head_u, double* head_t, char type) {
	if (type == 'C')   // f(1) - f(-1) / dt; head on -1
		return ((head_u[2] - head_u[0]) / (head_t[2] - head_t[0]));
	else if (type == 'F')  //-3f(0) + 4f(1) - f(2); head on 0
		return (-3 * head_u[0] + 4 * head_u[1] - head_u[2])
				/ (head_t[2] - head_t[0]);
	else
		// 3f(0) - 4f(-1) + f(-2)   head on -2
		return (3 * head_u[2] - 4 * head_u[1] + head_u[0])
				/ (head_t[2] - head_t[0]);
}

//*********************
//get initial pos of possible AT
//**********************
void Simulator::iniAT(int& s_pos, int& f_pos, int& rs_pos, int& rf_pos,
		int& max_pos, double& max, double* u) {
	int peak_pos, re_pos;
	finPeak(peak_pos, re_pos, u);

	if (peak_pos != step_save) {
		if (peak_pos - INTERVAL_AT < 0)
			s_pos = 0;
		else
			s_pos = peak_pos - INTERVAL_AT;
		if (peak_pos + INTERVAL_AT > step_save - 1)
			f_pos = step_save - 1;
		else
			f_pos = peak_pos + INTERVAL_AT;

		if (re_pos - INTERVAL_RT_F < 0)
			rs_pos = 0;
		else
			rs_pos = re_pos - INTERVAL_RT_F;
		if (re_pos + INTERVAL_RT_B > step_save - 1)
			rf_pos = step_save - 1;
		else
			rf_pos = re_pos + INTERVAL_RT_B;

		max_pos = s_pos;
		if (t)
			max = t[max_pos];
		else
			max = s_pos * dt_save;
	} else
		max_pos = step_save;
}

//*********************************
//find the peak of TMP wave
//***********************************
void Simulator::finPeak(int& peak_pos, int& re_pos, double* u) {
	peak_pos = 0;
	//time index of all plateau
	int* plateau = new int[step_save];
	for (int i = 0; i < step_save; i++)
		plateau[i] = 0;

	int pos = 0;
	for (int i = 0; i < step_save; i++) {
		if (u[i] >= THRESH_UP)
			plateau[pos++] = i;
	}

	if (pos == 0)         // no active 
			{
		peak_pos = step_save;
		re_pos = step_save;
	} else {
		// number of suspective plateaus
		int num_p = 1;
		for (int i = 1; i < pos; i++) {
			if ((plateau[i] - plateau[i - 1]) > 1)
				num_p++;
		}

		if (num_p == 1) {
			peak_pos = plateau[0];
			re_pos = plateau[pos - 1];
		} else {
			//the beginning idx for each plateau
			int* start_p = new int[num_p];
			int* len_p = new int[num_p];
			for (int i = 0; i < num_p; i++)
				len_p[i] = 0;

			start_p[0] = plateau[0];

			int j = 1;
			for (int i = 1; i < pos; i++) {
				if ((plateau[i] - plateau[i - 1]) == 1)
					len_p[j - 1]++;
				else
					start_p[j++] = plateau[i];
			}

			//pick the longest plateau
			int max_p = len_p[0];
			peak_pos = start_p[0];
			re_pos = peak_pos + len_p[0];
			for (int i = 1; i < num_p; i++) {
				if (len_p[i] > max_p) {
					max_p = len_p[i];
					peak_pos = start_p[i];
					re_pos = peak_pos + len_p[i];
				}
			}

			delete[] start_p;
			delete[] len_p;
		}
	}
	delete[] plateau;
}

//****************************************
//get AT & AT_seq using above functions
//****************************************
void Simulator::u2t() {
	double* div = new double[ID * step_save];
	double* u = (~state).getDataPtr();

	for (int i = 0; i < ID * step_save; i++)
		div[i] = 0;

	int pos = 0;
	double max;
	//numerical derivative
	for (int i = 0; i < ID; i++) {
		pos = 0;
		for (int j = 0; j < step_save; j++) {
			if (j == 0)
				div[pos] = numDiff(&u[pos], &t[j], 'F');
			else if (j == step_save - 1)
				div[pos] = numDiff(&u[pos - 2], &t[j - 2], 'B');
			else
				div[pos] = numDiff(&u[pos - 1], &t[j - 1], 'C');

			pos++;
		}
	}

	//find the maximum positiv
	int max_pos, s_pos, f_pos, rs_pos, rf_pos, j;
	for (int i = 0; i < ID; i++) {
		pos = i * step_save;
		iniAT(s_pos, f_pos, rs_pos, rf_pos, max_pos, max, &u[pos]);
		//determine AT ***********
		//not active nodes
		if (max_pos == step_save) {
			AT[i] = tmax;
			AT_idx[i] = step_save - 1;
			RT[i] = tmax;
			RT_idx[i] = step_save - 1;
		} else {
			for (j = s_pos; j <= f_pos; j++) {
				if (div[pos + j] > max) {
					max = div[pos + j];
					max_pos = j;
				}
			}
			AT[i] = t[max_pos];
			AT_idx[i] = max_pos;
			//determine repo time*******
			for (j = rs_pos; j < rf_pos; j++) {
				if (u[pos + j] <= u[pos + max_pos]) {
					RT[i] = t[j];
					RT_idx[i] = j;
					break;
				}
			}
			if (j == rf_pos) {
				RT[i] = tmax;
				RT_idx[i] = step_save - 1;
			}
		}
		APD[i] = RT[i] - AT[i];
	}

	delete[] div;
}

//********************************************************
//write output
//********************************************************
void Simulator::write(char* path) {
	char sid[5];
	sprintf(sid, "%d", ID);
	char* out = new char[MAX_FILENAME];
	char* file = new char[MAX_FILENAME];
	strcpy(out, path);
	strcat(out, "/Simulation/");
	strcat(out, sid);
	strcat(out, "/Output/");

	char subfolder[MAX_FILENAME];
	printf("Enter the subfolder for Simulation output path: \n");
	scanf("%s", &subfolder);
	strcat(out, subfolder);

	strcpy(file, out);
	strcat(file, "TMP.bin");
	state.save(file, "wb");

	strcpy(file, out);
	strcat(file, "BSP.bin");
	mea.save(file, "wb");

	strcpy(file, out);
	strcat(file, "AT.bin");
	writeFile(file, AT, sizeof(double), ID, "wb");

	strcpy(file, out);
	strcat(file, "AT.seq");
	writeFile(file, AT_idx, sizeof(int), ID, "wb");

	strcpy(file, path);
	strcat(file, "RT.bin");
	writeFile(file, RT, sizeof(double), ID, "wb");

	strcpy(file, out);
	strcat(file, "RT.seq");
	writeFile(file, RT_idx, sizeof(int), ID, "wb");

	strcpy(file, out);
	strcat(file, "APD.bin");
	writeFile(file, APD, sizeof(double), ID, "wb");

	strcpy(file, out);
	strcat(file, "Delay.sti");
	Sti_all.save(file, "wb");

	delete[] file;
	delete[] out;
}

