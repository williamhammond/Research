//MFree.cpp

//**************************************************************************************************************
//Inclusions
//**************************************************************************************************************
#include "inclu.h"
#include"BEMCell.h"
#include "external.h"
#include "MFree.h"
#include "BEM.h"
#include "Matrix.h"

//**************************************************************************************************************
//Definition of classes
//class: bgMesh
//**************************************************************************************************************
//*******************************************************
//ctors
//set the domain of the mesh as the first mesh , and set natural gp info for later use
// set the initial bgmesh as outside the original domain(i.e. isOut = true)
//*******************************************************
bgMesh::bgMesh(MFree& mfree) {
	xl = mfree.m_spBouCor[0].x;
	xu = xl + mfree.m_nMesInt;
	yl = mfree.m_spBouCor[0].y;
	yu = yl + mfree.m_nMesInt;
	zl = mfree.m_spBouCor[0].z;
	zu = zl + mfree.m_nMesInt;

	m_dpNatGPCor = new double[MAX_NAT_GP];
	m_dpNatGPWei = new double[MAX_NAT_GP];

	m_nSta = 2;
	m_nGPNum = 0;
}

//*******************************************************
//dtor
//*******************************************************
bgMesh::~bgMesh() {

	delete[] m_dpNatGPCor;
	delete[] m_dpNatGPWei;
}

//*******************************************************
//setMesSeq
//change the mesh domain as loop continues, i.e. to sweep through all the meshes 
//*******************************************************
void bgMesh::setMesSeq(int i, int j, int k, MFree& mfree) {
	// generate domain of specified mesh
	xl = mfree.m_spBouCor[0].x + i * mfree.m_nMesInt;
	xu = xl + mfree.m_nMesInt;

	yl = mfree.m_spBouCor[0].y + j * mfree.m_nMesInt;
	yu = yl + mfree.m_nMesInt;

	zl = mfree.m_spBouCor[0].z + k * mfree.m_nMesInt;
	zu = zl + mfree.m_nMesInt;

	//whether it's out
	comSta(mfree);
}

//*******************************************************
//isout
//determine whether the mesh is totally inside the original domain
//set the "m_nSta", with the whole coordinate info from MFree class
// m_nSta = 0, while totally out
// =1, when partially out, in such case, gps' position needs to be checked
// else, totally in 
//*******************************************************
void bgMesh::comSta(MFree& mfree) {
	m_nSta = 2;
	m_nGPNum = 0; // re-initialize

	for (int i = 0; i < mfree.m_nNodNum; i++) {
		//at least one nodal is in the mesh, so it's not totally out
		if (mfree.m_spNodCor[i].x >= xl && mfree.m_spNodCor[i].x <= xu
				&& mfree.m_spNodCor[i].y >= yl && mfree.m_spNodCor[i].y <= yu
				&& mfree.m_spNodCor[i].z >= zl && mfree.m_spNodCor[i].z <= zu)

			m_nGPNum = m_nGPNum + 1;
	}

	if (m_nGPNum == 0)
		m_nSta = 0;

	// generate natural gp coordinates&weights, just for 1D info

	else {
		m_nGPNum = floor(pow(m_nGPNum, 1.0 / 3.0) + 2);
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
			break;
		}
	}

}

//*******************************************************
//genGloGP
//generate global gp coordinates and weights, one point each time,
//loop through all the gps in a mesh
//*******************************************************
void bgMesh::genGloGP(cor &gpcor, cor& gpwei, int i, int j, int k) {
	// x = ( xl + xu) /2 + gp*(xu-xl)/2
	gpcor.x = (xl + xu) / 2.0 + m_dpNatGPCor[i] * (xu - xl) / 2.0;
	gpcor.y = (yl + yu) / 2.0 + m_dpNatGPCor[j] * (yu - yl) / 2.0;
	gpcor.z = (zl + zu) / 2.0 + m_dpNatGPCor[k] * (zu - zl) / 2.0;
	// neww = (xu-xl)*oldw/2
	gpwei.x = m_dpNatGPWei[i] * (xu - xl) / 2.0;
	gpwei.y = m_dpNatGPWei[j] * (yu - yl) / 2.0;
	gpwei.z = m_dpNatGPWei[k] * (zu - zl) / 2.0;
}

//**************************************************************************************************************-
//class:nodal
//**************************************************************************************************************
//*******************************************************
//ctors
//get the cor of the nodal according to input seq, and generate p for computing shape function
//*******************************************************
nodal::nodal() :
		p(4, 1) {
//	p = new double[4];
	// no assignment, just generate the structure
}

//*******************************************************
//dtor
//*******************************************************
nodal::~nodal() {
//	delete p;
}

//*******************************************************
//reIni
//get nodal info for each nodal according to its sequence
//*******************************************************-------
void nodal::reIni(int seq, MFree& mfree) {
	m_cor = mfree.m_spNodCor[seq - 1];

//	p = vertcat(1, m_cor.x, m_cor.y, m_cor.z,m_cor.x*m_cor.y, m_cor.x*m_cor.z, m_cor.y*m_cor.z,
//	m_cor.x*m_cor.x, m_cor.y*m_cor.y, m_cor.z*m_cor.z);
	p.setData(1, 0, 0);
	p.setData(m_cor.x, 1, 0);
	p.setData(m_cor.y, 2, 0);
	p.setData(m_cor.z, 3, 0);
}

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

MatrixGroup::~MatrixGroup() {
}

//**************************************************************************************************************
//class: gp --- Important class
//**************************************************************************************************************
//*******************************************************
//ctors
//only pre-allocate storage, the specific info is got through re-initialize function
// for every gp
//*******************************************************
//without information of number of influence nodes************************************
//sha & B are empty matrix, and wei are empty*****************************************
//useful for single gp operations******************************************************
gp::gp() :
		D0(3, 3), D(3, 3), Q(3, 3), shape(), p(4, 1), sha(), A(4, 4), B(), invA(
				4, 4) {
	Q.setData(1.0, 0, 0);
	Q.setData(1.0, 1, 1);
	Q.setData(1.0, 2, 2);
	if (ISO == false) {
		D0.setData(DI_S, 0, 0);
		D0.setData(DI_S, 1, 1);
		D0.setData(DI_L, 2, 2);
	} else {
		D0.setData(DI, 0, 0);
		D0.setData(DI, 1, 1);
		D0.setData(DI, 2, 2);
	}

	m_ptrNod = new nodal();

	m_nSur = 0;
	m_npSurSeq = NULL;
	m_nodWei = NULL;
	m_nodWeiDx = NULL;
	m_nodWeiDy = NULL;
	m_nodWeiDz = NULL;

	m_nCurrent = 0;
}

//****************************************************
//for a bunch of gpts, use fixed length to initialize all memory related to 
//influence nodes (use MAX_INF)
//****************************************************
gp::gp(int num) :
		D0(3, 3), D(3, 3), Q(3, 3), shape(3, num), p(4, 1), sha(1, num), A(4,
				4), B(4, num), invA(4, 4) {
	Q.setData(1.0, 0, 0);
	Q.setData(1.0, 1, 1);
	Q.setData(1.0, 2, 2);
	if (ISO == false) {
		D0.setData(DI_S, 0, 0);
		D0.setData(DI_S, 1, 1);
		D0.setData(DI_L, 2, 2);
	} else {
		D0.setData(DI, 0, 0);
		D0.setData(DI, 1, 1);
		D0.setData(DI, 2, 2);
	}

	m_ptrNod = new nodal();

	m_nSur = 0;
	m_npSurSeq = new int[num];
	m_nodWei = new double[num];
	m_nodWeiDx = new double[num];
	m_nodWeiDy = new double[num];
	m_nodWeiDz = new double[num];

	m_nCurrent = 0;
}

//*******************************************************
//dtor
//*******************************************************
gp::~gp() {
	delete[] m_npSurSeq;
	delete[] m_nodWei;
	delete[] m_nodWeiDx;
	delete[] m_nodWeiDy;
	delete[] m_nodWeiDz;

	delete m_ptrNod;
}

//*******************************************************
//reIni
//re-initialize, i.e. re-set the info of gp when loop through each gp
//*******************************************************
void gp::reIni(cor& gpCor, cor& gpWei, MFree& mfree) {
	m_gpCor = gpCor;
	m_dWei = gpWei.x * gpWei.y * gpWei.z;
	m_gpFib.x = m_gpFib.y = m_gpFib.z = 0;
	m_nCurrent++;

	Poly3D(p.d, m_gpCor.x, m_gpCor.y, m_gpCor.z);
	Poly3D_dx(p.dx, m_gpCor.x, m_gpCor.y, m_gpCor.z);
	Poly3D_dy(p.dy, m_gpCor.x, m_gpCor.y, m_gpCor.z);
	Poly3D_dz(p.dz, m_gpCor.x, m_gpCor.y, m_gpCor.z);

	sha.reSet();
	A.reSet();
	B.reSet();
	invA.reSet();
	shape.reSet(0.0);

	getSur(mfree);

}

//*******************************************************
//getSur
// get the supportive nodals for each gp, record its # and sequence
// this function is called when each new gp is dealt with, and no need to 
// be applied within the same gp when dealing with different supportive nodals
//*******************************************************
void gp::getSur(MFree& mfree) {
	double a, r, dx, dy, dz;
	m_nSur = 0;
	for (int k = 0; k < MAX_INF; k++) {
		m_npSurSeq[k] = 0;
		m_nodWei[k] = 0;
		m_nodWeiDx[k] = 0;
		m_nodWeiDy[k] = 0;
		m_nodWeiDz[k] = 0;
	}
	//loop through each nodal, to check which ones influence domain cover gp
	for (int i = 0; i < mfree.m_nNodNum; i++) { //the distance between gp and current nodals
		a = DMAX * mfree.m_dpNodDom[i];
		dx = mfree.m_spNodCor[i].x - m_gpCor.x;
		dy = mfree.m_spNodCor[i].y - m_gpCor.y;
		dz = mfree.m_spNodCor[i].z - m_gpCor.z;
		r = sqrt(dx * dx + dy * dy + dz * dz);

		if (r > 0.00001) {
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
					m_nodWei[m_nSur] = 2.0 / 3.0 - 4 * r * r + 4 * r * r * r;
					// wx = -8rx + 12r*rx    ( rx = (x - xi)/(a^2*r))
					m_nodWeiDx[m_nSur] = (-8 * r + 12 * r * r) * dx;
					m_nodWeiDy[m_nSur] = (-8 * r + 12 * r * r) * dy;
					m_nodWeiDz[m_nSur] = (-8 * r + 12 * r * r) * dz;
					// wxx = -8 + 12(rx^2/r + r)
//		m_nodWeiDD.x = (-8+24*r)*dx*dx + m_nodWeiD.x*ddx/dx;
//		m_nodWeiDD.y = (-8+24*r)*dy*dy + m_nodWeiD.x*ddy/dy;
//		m_nodWeiDD.z = (-8+24*r)*dz*dz + m_nodWeiD.z*ddz/dz;
					m_npSurSeq[m_nSur++] = i + 1;
				}
			} else if (0.5 < r && r <= 1) {
				if (m_nSur + 1 > MAX_INF) {
					cout << "influencing nodes out of limits" << endl;
					wait();
					exit(0);
				} else {
					// w = 4/3 - 4r + 4r^2 - 4/3r^3
					m_nodWei[m_nSur] = 4.0 / 3.0 - 4 * r + 4 * r * r
							- 4.0 / 3.0 * r * r * r;
					// wx = -4rx/r + 8rx - 4rrx
					m_nodWeiDx[m_nSur] = (-4 + 8 * r - 4 * r * r) * dx;
					m_nodWeiDy[m_nSur] = (-4 + 8 * r - 4 * r * r) * dy;
					m_nodWeiDz[m_nSur] = (-4 + 8 * r - 4 * r * r) * dz;
					//wxx = 8 -4(1/r - rx^2/r^3) - 4(rx^2 + r)
//		m_nodWeiDD.x = (8 -8*r)*dx*dx + m_nodWeiD.x*ddx/dx;
//		m_nodWeiDD.y = (8 -8*r)*dy*dy + m_nodWeiD.y*ddy/dy;
//		m_nodWeiDD.z = (8 -8*r)*dz*dz + m_nodWeiD.z*ddz/dz;
					m_npSurSeq[m_nSur++] = i + 1;
				}
			}
		}
	}
//	cout<<"Number of support nodes: "<<m_nSur<<endl;
}

//*******************************************************
//sinTrans
//for each gp, calculate its single transfer matrix, i.e. f(xi), by adding A,B for every ( gp,nodal) pair
//according to "matrix calculus"
// if flag == 1, new method, i.e. assemble shape function only once, and deal with all boundary points
//*******************************************************-  
void gp::sinTrans(MFree& mfree) {
//	cout<<"Generating matrix A & B"<<endl;

	Matrix tmp(4, 4);

	//loop through every supportive point
	//construct A,B, and their derivatives
	for (int i = 0; i < m_nSur; i++) {
		//compute A,B and their derivatives nodal by nodal
		m_ptrNod->reIni(m_npSurSeq[i], mfree);
		tmp = m_ptrNod->p * (~(m_ptrNod->p));
		A.d += tmp * m_nodWei[i];
		A.dx += tmp * m_nodWeiDx[i];
		A.dy += tmp * m_nodWeiDy[i];
		A.dz += tmp * m_nodWeiDz[i];

		B.d.setSub((m_ptrNod->p) * m_nodWei[i], 0, 3, i, i);
		B.dx.setSub((m_ptrNod->p) * m_nodWeiDx[i], 0, 3, i, i);
		B.dy.setSub((m_ptrNod->p) * m_nodWeiDy[i], 0, 3, i, i);
		B.dz.setSub((m_ptrNod->p) * m_nodWeiDz[i], 0, 3, i, i);
	}
	//compute ddSha

	invA.d = A.d.inv();

//	cout<<"Generating shape function"<<endl;

	sha.d.setSub((~(p.d)) * invA.d * B.d.sub(0, 3, 0, m_nSur - 1), 0, 0, 0,
			m_nSur - 1);
//	sha.d.save("../BEMForward/GP/sha.bin","a+b");
	if (checkNan(sha.d)) {
		cout << "Nan happens @ sha" << endl;
		wait();
	}

	int dsum = sha.d.sum() - 1.0;
	if (dsum > 0.00001 || dsum < -0.00001) {
		cout << "Error: sum of shape is not equal to 1 on gpt" << endl;
		wait();
	}

	invA.dx = -invA.d * A.dx * invA.d;
//	invAdd[0] = -(invAd[0]*Ad[0]*invA + invA*Add[0]*invA + invA*Ad[0]*invAd[0]);
	invA.dy = -invA.d * A.dy * invA.d;
//	invAdd[1] = -(invAd[1]*Ad[1]*invA + invA*Add[1]*invA + invA*Ad[1]*invAd[1]);
	invA.dz = -invA.d * A.dz * invA.d;
//	invAdd[2] = -(invAd[2]*Ad[2]*invA + invA*Add[2]*invA + invA*Ad[2]*invAd[2]);

	/*	if (checkNan(invA.d))
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
	sha.dx.setSub(
			(~p.dx) * invA.d * B.d.sub(0, 3, 0, m_nSur - 1)
					+ (~p.d) * invA.dx * B.d.sub(0, 3, 0, m_nSur - 1)
					+ (~p.d) * invA.d * B.dx.sub(0, 3, 0, m_nSur - 1), 0, 0, 0,
			m_nSur - 1);
	sha.dy.setSub(
			(~p.dy) * invA.d * B.d.sub(0, 3, 0, m_nSur - 1)
					+ (~p.d) * invA.dy * B.d.sub(0, 3, 0, m_nSur - 1)
					+ (~p.d) * invA.d * B.dy.sub(0, 3, 0, m_nSur - 1), 0, 0, 0,
			m_nSur - 1);
	sha.dz.setSub(
			(~p.dz) * invA.d * B.d.sub(0, 3, 0, m_nSur - 1)
					+ (~p.d) * invA.dz * B.d.sub(0, 3, 0, m_nSur - 1)
					+ (~p.d) * invA.d * B.dz.sub(0, 3, 0, m_nSur - 1), 0, 0, 0,
			m_nSur - 1);

	if (checkNan(sha.dx)) {
		cout << "Nan happens @ sha_dx" << endl;
		wait();
	}
	if (checkNan(sha.dy)) {
		cout << "Nan happens @ sha_dy" << endl;
		wait();
	}
	if (checkNan(sha.dz)) {
		cout << "Nan happens @ sha_dz" << endl;
		wait();
	}
	/*	ddSha[0] = 2*dp[0]*(invAd[0]*B + invA*Bd[0]) + p*(invAdd[0]*B + invA*Bdd[0] + 2*invAd[0]*Bd[0]) + ddp[0]*invA*B;
	 ddSha[1] = 2*dp[1]*(invAd[1]*B + invA*Bd[1]) + p*(invAdd[1]*B + invA*Bdd[1] + 2*invAd[1]*Bd[1])+ ddp[1]*invA*B;
	 ddSha[2] = 2*dp[2]*(invAd[2]*B + invA*Bd[2]) + p*(invAdd[2]*B + invA*Bdd[2] + 2*invAd[2]*Bd[2])+ ddp[2]*invA*B;
	 */
	/*	mwArray shape = vertcat(ddSha[0], ddSha[1], ddSha[2]);



	 //here to interpolate gp's con ---- using similar idea of mfree approximation
	 m_sinTrans = m_dpR*m_dpCon*shape;
	 */
//	cout<<"Generating fiber for gs"<<endl;
	if (ISO == 0)
		setGPCon(mfree, NULL);

//	cout<<"Set shape functions"<<endl;

	shape.setSub(sha.dx.sub(0, 0, 0, m_nSur - 1), 0, 0, 0, m_nSur - 1);
	shape.setSub(sha.dy.sub(0, 0, 0, m_nSur - 1), 1, 1, 0, m_nSur - 1);
	shape.setSub(sha.dz.sub(0, 0, 0, m_nSur - 1), 2, 2, 0, m_nSur - 1);

//	mwArray shape = vertcat(ddSha[0], ddSha[1], ddSha[2]);
//	double* ptr = shape.getDataPtr();
//	writeFile("../BEMForward/GP/shapeD.bin",ptr,sizeof(double),3*m_nSur,"a+b");

	//***store into mfree class************
//	cout<<"Save gpts"<<endl;

	if (mfree.m_nGPNum + 1 > MAX_GP) {
		cout << "not enough memeory for gps: " << endl;
		exit(0);

	}

	//cout<<"Store gpts information"<<endl;

	mfree.m_spGPCor[mfree.m_nGPNum] = m_gpCor; // gp cor,record after first sweep
	mfree.m_spGPFib[mfree.m_nGPNum] = m_gpFib;
	mfree.m_dpGPWei[mfree.m_nGPNum] = m_dWei;
	mfree.m_npGPNBNum[mfree.m_nGPNum] = m_nSur; //number of NB nodes for each gp
	mfree.m_mpGPshaD[mfree.m_nGPNum].setData(0.0, 3, m_nSur);
	mfree.m_mpGPshaD[mfree.m_nGPNum] = D * (shape.sub(0, 2, 0, m_nSur - 1))
			* m_dWei;

	mfree.m_mpGPsha[mfree.m_nGPNum].setData(0.0, 4, m_nSur);
	mfree.m_mpGPsha[mfree.m_nGPNum].setSub(sha.d.sub(0, 0, 0, m_nSur - 1), 0, 0,
			0, m_nSur - 1);
	mfree.m_mpGPsha[mfree.m_nGPNum].setSub(shape.sub(0, 2, 0, m_nSur - 1), 1, 3,
			0, m_nSur - 1);

	mfree.m_nGPNum++;

	for (int j = 0; j < m_nSur; j++) {
		mfree.m_npGPNB[mfree.m_nCurrentGP] = m_npSurSeq[j]; //sequnce of NB nodes for each gp
//		mfree.m_dpGPsha[mfree.m_nCurrentGP] = sha.d(0,j);
		mfree.m_nCurrentGP++;
	}

}

//*******************************************************
//sin2Glo
//assemble single transfer vector( only includes non-zeros) into global transfer vector( including all nodals)
//input: gloTrans is the current row (corresponding to current surface points) and can be modified by this function
//*******************************************************
/*void gp::sin2Glo(mwArray& gloTrans, int nodSeq )
 {
 for(int i=0; i<m_nSur; i++)
 gloTrans(nodSeq,m_npSurSeq[i]) =gloTrans(nodSeq,m_npSurSeq[i])+ m_sinTrans(i+1);
 }
 */

//***************************************************
//computeQ: compute the transform matrix between two coordinate system
// global = Q local      (on basis sense)
// Qij: i  for local, j for global
//***************************************************
void gp::computeQ() {
	double e1, e2, e3;
	if (m_gpFib.z != 0 && m_gpFib.x != 0) {
		//local x
		e1 = m_gpFib.x;
		e2 = m_gpFib.y;
		e3 = -(m_gpFib.x * m_gpFib.x + m_gpFib.y * m_gpFib.y) / (m_gpFib.z);
		if (m_gpFib.z < 0) {
			e1 = -e1;
			e2 = -e2;
			e3 = -e3;
		}
		Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
		Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
		Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
		//for local y axis
		e1 = -m_gpFib.y * m_gpFib.y / m_gpFib.x;
		e2 = m_gpFib.y;
		e3 = 0;
		if (m_gpFib.y * m_gpFib.x < 0) {
			e1 = -e1;
			e2 = -e2;
			e3 = -e3;
		}
		Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
		Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
		Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
		//for local z axis (longitudal fiber)
		e1 = m_gpFib.x;
		e2 = m_gpFib.y;
		e3 = m_gpFib.z;
		Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
		Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
		Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 2);
	} else {
		//fiber direction in x-y plane
		if (m_gpFib.z == 0 && m_gpFib.x != 0) {
			e1 = 0;
			e2 = 0;
			e3 = -1;
			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
			//for local y axis
			e1 = -m_gpFib.y * m_gpFib.y / m_gpFib.x;
			e2 = m_gpFib.y;
			e3 = 0;
			if (m_gpFib.x * m_gpFib.y < 0) {
				e1 = -e1;
				e2 = -e2;
				e3 = -e3;
			}
			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
			//for local z axis (longitudal fiber)
			e1 = m_gpFib.x;
			e2 = m_gpFib.y;
			e3 = m_gpFib.z;
			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 2);
		}
		//fiber direction in y-z plane
		else if (m_gpFib.x == 0 && m_gpFib.z != 0) {
			e1 = m_gpFib.x;
			e2 = m_gpFib.y;
			e3 = -(m_gpFib.x * m_gpFib.x + m_gpFib.y * m_gpFib.y) / (m_gpFib.z);
			if (m_gpFib.y * m_gpFib.z < 0) {
				e1 = -e1;
				e2 = -e2;
				e3 = -e3;
			}
			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
			//for local y axis
			e1 = -1;
			e2 = 0;
			e3 = 0;

			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
			//for local z axis (longitudal fiber)
			e1 = m_gpFib.x;
			e2 = m_gpFib.y;
			e3 = m_gpFib.z;
			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 1);
		} else {
			if (m_gpFib.y > 0) {
				e1 = 0;
				e2 = 0;
				e3 = -1;
			} else {
				e1 = 0;
				e2 = 0;
				e3 = 1;
			}
			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 0, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 0, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 0, 2);
			//for local y axis
			e1 = -1;
			e2 = 0;
			e3 = 0;

			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 1, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 1, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 1, 2);
			//for local z axis (longitudal fiber)
			e1 = m_gpFib.x;
			e2 = m_gpFib.y;
			e3 = m_gpFib.z;
			Q.setData(::dotProduct(1, 0, 0, e1, e2, e3), 2, 0);
			Q.setData(::dotProduct(0, 1, 0, e1, e2, e3), 2, 1);
			Q.setData(::dotProduct(0, 0, 1, e1, e2, e3), 2, 2);
		}
	}
}

//*******************************************************
//setGPCon
//get the conductivity of gps using interpolation  --- 3*3 tensor
//or closest points to reference heart geometry
//or directly read in if existed
//*******************************************************
void gp::setGPCon(MFree& mfree, double* con) {

	//interpolated from mfree particles
	if (FIBER_INTER) {
		Matrix fSur_x(m_nSur, 1);
		Matrix fSur_y(m_nSur, 1);
		Matrix fSur_z(m_nSur, 1);
		for (int i = 0; i < m_nSur; i++)
			fSur_x.setData(mfree.m_spNodCon[m_npSurSeq[i] - 1].x, i, 0);
		for (int i = 0; i < m_nSur; i++)
			fSur_y.setData(mfree.m_spNodCon[m_npSurSeq[i] - 1].y, i, 0);
		for (int i = 0; i < m_nSur; i++)
			fSur_z.setData(mfree.m_spNodCon[m_npSurSeq[i] - 1].z, i, 0);

		m_gpFib.x = (sha.d.sub(0, 0, 0, m_nSur - 1) * fSur_x)(0, 0);
		m_gpFib.y = (sha.d.sub(0, 0, 0, m_nSur - 1) * fSur_y)(0, 0);
		m_gpFib.z = (sha.d.sub(0, 0, 0, m_nSur - 1) * fSur_z)(0, 0);
	} else if (FIBER_CLOSE) {
		double* ref_cor = NULL, *ref_fib = NULL;
		int ref_num = 0;
		readFile(FILE_REF_COR, &ref_cor, sizeof(double), 3, ref_num);
		readFile(FILE_REF_FIB, &ref_fib, sizeof(double), 3, ref_num);

		double dx, dy, dz, r, r_min;
		int i_min = 0, count = 0;

		dx = m_gpCor.x - ref_cor[count++];
		dy = m_gpCor.y - ref_cor[count++];
		dz = m_gpCor.z - ref_cor[count++];
		r_min = dx * dx + dy * dy + dz * dz;

		for (int i = 1; i < ref_num; i++) {
			dx = m_gpCor.x - ref_cor[count++];
			dy = m_gpCor.y - ref_cor[count++];
			dz = m_gpCor.z - ref_cor[count++];

			r = dx * dx + dy * dy + dz * dz;

			if (r_min > r) {
				r_min = r;
				i_min = i;
			}
		}

		count = i_min * DIM;
		m_gpFib.x = ref_fib[count];
		m_gpFib.y = ref_fib[count + 1];
		m_gpFib.z = ref_fib[count = 2];

	} else {
		m_gpFib.x = mfree.m_spGPFib[m_nCurrent - 1].x;
		m_gpFib.y = mfree.m_spGPFib[m_nCurrent - 1].y;
		m_gpFib.z = mfree.m_spGPFib[m_nCurrent - 1].z;
	}

	computeQ();
//	Q.save("../BEMForward/Test/Q.bin","wb");

	D = (~Q) * D0 * Q;
//	m_dpCon = transpose(Q)*D0*Q;
	/*	if(con)
	 {
	 mwArray conGP = zeros(1,2);
	 mwArray conSur = zeros(2,m_nSur);
	 pos = 1;
	 for ( int i=0; i< m_nSur; i++)
	 {
	 conSur(pos++) = con[m_npSurSeq[i]-1];
	 conSur(pos++) = con[m_npSurSeq[i]-1+mfree.m_nNodNum];
	 }
	 conGP = sha*transpose(conSur);
	 D0(1,1) = conGP.ExtractScalar(1);
	 D0(2,2) = D0(1,1);
	 D0(3,3) = conGP.ExtractScalar(2);
	 }
	 */
	//***********wrong****************************
	/* 	double rot_y = PI/2;
	 if(fx!=0)
	 rot_y = atan(fx / fy);

	 double rot_z = atan( sqrt(fx*fx+fy*fy)/fz);

	 double dRxy[] = {cos(rot_y), -sin(rot_y), 0, sin(rot_y), cos(rot_y), 0, 0, 0, 1};
	 mwArray  Rxy(3,3,dRxy);
	 double dRxz[] = {cos(rot_z), 0, -sin(rot_z), 0 ,1 ,0 , sin(rot_z), 0, cos(rot_z)};
	 mwArray Rxz(3,3,dRxz);

	 m_dpCon =  inv(Rxz*Rxy)*D0*Rxz*Rxy;
	 *//**/
	//**************read from gp.fib, real fiber directions of gps

//	delete []dp_fSur;
}

//**************************************************************************************************************
//class: MFree --- class providing interface to external application
//Read in data, which is further used by embeded classes
//Output data: 1, transfer matrix 2, solution given boundary values
//implement the overall "meshfree" step 
//**************************************************************************************************************
//*******************************************************
//ctors
// initialize: input # gps for each mesh and mesh interval, and use to generate bgmesh and assign the value to 
//             corresponding variables 
//             read in data about nodal info(cor,con)
//             compute the influence domain of every nodal
//*******************************************************
MFree::MFree(double mesInt, char** filNam, BEM& bem) :
		m_nMesInt(mesInt) {
	//read in data-------------------------------------------------------------
	//nodal cor

	//open file, and get the size of the file, then correpondingly # of data and further # of nodals
	//nodal cor
	double* data = NULL;
	int i, pos = 0;
	readFile(filNam[0], &data, sizeof(double), 3, m_nNodNum);
	m_spNodCor = new cor[m_nNodNum];
	//read in nodal by nodal
	for (i = 0; i < m_nNodNum; i++) {
		m_spNodCor[i].x = data[pos++];
		m_spNodCor[i].y = data[pos++];
		m_spNodCor[i].z = data[pos++];
	}
	if (data)
		delete[] data;
	data = NULL;
	//nodal con
	if (ISO == false) {
		readFile(filNam[1], &data, sizeof(double), 3, m_nNodNum);
		m_spNodCon = new cor[m_nNodNum];
		//read in nodal by nodal
		pos = 0;
		for (i = 0; i < m_nNodNum; i++) {
			m_spNodCon[i].x = data[pos++];
			m_spNodCon[i].y = data[pos++];
			m_spNodCon[i].z = data[pos++];
		}
		if (data)
			delete[] data;
		data = NULL;
	}

	readFile(filNam[2], &data, sizeof(double), 3, m_nTetNodNum);
	m_spTetCor = new cor[m_nTetNodNum];

	pos = 0;
	//read in nodal by nodal
	for (i = 0; i < m_nTetNodNum; i++) {
		m_spTetCor[i].x = data[pos++];
		m_spTetCor[i].y = data[pos++];
		m_spTetCor[i].z = data[pos++];
	}
	if (data)
		delete[] data;
	data = NULL;

	m_npTetSeq = NULL;
	readFile(filNam[3], &m_npTetSeq, sizeof(int), 4, m_nTetNum);

	//*****************************
	//***GP info allocated later, dependent on whether generate new or not****************
	m_spGPCor = NULL;
	m_spGPFib = NULL;
	m_dpGPWei = NULL;
	m_nGPNum = 0;
	m_nCurrentGP = 0;
	m_npGPNB = NULL;      //sequnce of NB nodes for each gp
	m_npGPNBNum = NULL;  //number of NB nodes for each gp
	m_mpGPsha = NULL;
	m_mpGPshaD = NULL; //[MAX_GP*MAX_INF*3];
	m_ptrGP = NULL;
//	m_dpGPsha = new double[MAX_GP*MAX_INF];

	m_dpR = new double[3];
	m_dpTrans = new double[m_nNodNum * (bem.m_nNodNum + 1)]; // transpose matrix first, easy to assign
	m_dpTrans_arranged = new double[m_nNodNum
			* (bem.m_nNodNum + 1 - bem.m_nForNum)];

	//compute the info of influence domain-----------------------------
	getNodInfDom(INFNUM);

	m_npMesNum = new int[4];
	genBgMes();
	m_ptrMes = new bgMesh(*this);
	//generate pointer to gp, without any specific info,just allocate room
}

//*******************************************************
//dtor
//*******************************************************
MFree::~MFree() {
	delete[] m_spNodCor;
	if (m_spNodCon)
		delete[] m_spNodCon;
	delete[] m_dpNodDom;
	delete[] m_spTetCor;
	delete[] m_npTetSeq;

	delete[] m_spGPCor;
	if (m_spGPFib)
		delete[] m_spGPFib;
	if (m_dpGPWei)
		delete[] m_dpGPWei;
	delete[] m_npGPNB;
	delete[] m_npGPNBNum;
	delete[] m_mpGPsha;
	delete[] m_mpGPshaD;

	delete[] m_spBouCor;
	delete[] m_npMesNum;
	delete[] m_dpR;
	delete[] m_dpTrans;
	delete[] m_dpTrans_arranged;
	delete m_ptrMes;
}

//*******************************************************
//getNodInfDom:
//according to specified # of neighbours for each nodal
//*******************************************************
void MFree::getNodInfDom(int num) {
	m_dpNodDom = new double[m_nNodNum];
	double* inf = new double[num]; // record of influence neighbors in ascent sequence
	double r = 0;
	for (int l = 0; l < m_nNodNum; l++) {
		for (int i = 0; i < num; i++)
			inf[i] = MAX_R;
		//loop through every nodal
		for (int j = 0; j < m_nNodNum; j++) {
			//loop through every nodal
			if (l != j) {
				r = sqrt(
						power(m_spNodCor[l].x - m_spNodCor[j].x)
								+ power(m_spNodCor[l].y - m_spNodCor[j].y)
								+ power(m_spNodCor[l].z - m_spNodCor[j].z));
				//check through the record queue, and add/update accordingly
				//determine the pos in record queue
				int k = 0;
				for (k = 0; k < num && r > inf[k]; k++) {
				}
				//update record queue( insert in new neighbor at correct pos, and move following neighbors backward)
				if (k < num) {
					for (int m = (num - 1); m > k; m--)
						inf[m] = inf[m - 1];
					inf[k] = r;
				}
			}
		}
		//get radius for influence domain
		//m_dpNodDom[l] = inf[num-1];
		m_dpNodDom[l] = 0;
		for (int j = 0; j < num; j++)
			m_dpNodDom[l] += inf[j];
		m_dpNodDom[l] /= num;
	}
	delete[] inf;
}

//*******************************************************
//genBgMes
//generate bg mesh according to read-in data of nodals, actually target is to get the boundary cor for mesh, and
//basic info about mesh, eg. #gp, #meshes
//*******************************************************
void MFree::genBgMes() {
	//pre-allocate storage
	m_spBouCor = new cor[2];

	m_spBouCor[0] = m_spNodCor[0];
	m_spBouCor[1] = m_spNodCor[0];

	//loop through every tetra nodes to determine the lower bound and upper bound of bg mesh
	for (int i = 0; i < m_nTetNodNum; i++) {
		if (m_spBouCor[0].x > m_spTetCor[i].x)
			m_spBouCor[0].x = m_spTetCor[i].x;
		if (m_spBouCor[0].y > m_spTetCor[i].y)
			m_spBouCor[0].y = m_spTetCor[i].y;
		if (m_spBouCor[0].z > m_spTetCor[i].z)
			m_spBouCor[0].z = m_spTetCor[i].z;
		if (m_spBouCor[1].x < m_spTetCor[i].x)
			m_spBouCor[1].x = m_spTetCor[i].x;
		if (m_spBouCor[1].y < m_spTetCor[i].y)
			m_spBouCor[1].y = m_spTetCor[i].y;
		if (m_spBouCor[1].z < m_spTetCor[i].z)
			m_spBouCor[1].z = m_spTetCor[i].z;
	}
	//compute the total# mesh
	//get the directional #  -- round off into higher integral
	double dx = (m_spBouCor[1].x - m_spBouCor[0].x);
	double dy = (m_spBouCor[1].y - m_spBouCor[0].y);
	double dz = (m_spBouCor[1].z - m_spBouCor[0].z);
	m_npMesNum[1] = ceil(dx / m_nMesInt);
	m_npMesNum[2] = ceil(dy / m_nMesInt);
	m_npMesNum[3] = ceil(dz / m_nMesInt);
	// adjust the boundary cor to make integral# mesh
	dx = m_npMesNum[1] * m_nMesInt - dx;
	dy = m_npMesNum[2] * m_nMesInt - dy;
	dz = m_npMesNum[3] * m_nMesInt - dz;
	m_spBouCor[0].x -= dx / 2;
	m_spBouCor[0].y -= dy / 2;
	m_spBouCor[0].z -= dz / 2;
	m_spBouCor[1].x += dx / 2;
	m_spBouCor[1].y += dy / 2;
	m_spBouCor[1].z += dz / 2;
	//total #
	m_npMesNum[0] = m_npMesNum[1] * m_npMesNum[2] * m_npMesNum[3];
	cout << "Number of background meshes: " << m_npMesNum[0];
	//total # of gp
	//m_nGPTot = power(m_nGPNum)*m_nGPNum * m_npMesNum[0];
}

//*******************************************************
//setSurPoiInf
//set the current surface point info, which is got from external application, here BEM object
//*******************************************************
void MFree::setSurPoiInf(int seq, cor poiCor)	//, cor poiCon)
		{
	m_nSurPoiSeq = seq;
	m_sSurPoiCor = poiCor;
//	m_dSurPoiCon = poiCon.x;
}
//*******************************************************
//isout
//judge whether gp is out of volume
//*******************************************************
bool MFree::isout(cor gpCor) {
	//loop through all tetras, to find whether at least one include this gp
	int pos = 0;
	for (int i = 0; i < m_nTetNum; i++) {
		int n0 = m_npTetSeq[pos++];
		int n1 = m_npTetSeq[pos++];
		int n2 = m_npTetSeq[pos++];
		int n3 = m_npTetSeq[pos++];

		cor c0 = m_spTetCor[n0 - 1];
		cor c1 = m_spTetCor[n1 - 1];
		cor c2 = m_spTetCor[n2 - 1];
		cor c3 = m_spTetCor[n3 - 1];

		int V0 = inVol(c0, c1, c2, c3);
		int V1 = inVol(gpCor, c1, c2, c3);
		int V2 = inVol(c0, gpCor, c2, c3);
		int V3 = inVol(c0, c1, gpCor, c3);
		int V4 = inVol(c0, c1, c2, gpCor);

		if (V0 == 0) {
			//	cerr << "The tetrahedron is degenerate !!" << endl;
			//	exit(0);
			//	return true;
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

//***************************************************************************
//writeGPInfo: record all in gp info, including cor / fib / inf_num / inf_seq /wei /shaD
//all gpt info are saved in //GP// directory
//**************************************************************************
void MFree::writeGPInfo(char* path) {
	char* file = new char[MAX_FILENAME];
	int count = 0;

	strcpy(file, path);
	strcat(file, "gp.cor");
	double* gp = new double[3 * m_nGPNum];
	for (int j = 0; j < m_nGPNum; j++) {
		gp[count++] = m_spGPCor[j].x;
		gp[count++] = m_spGPCor[j].y;
		gp[count++] = m_spGPCor[j].z;
	}
	writeFile(file, gp, sizeof(double), 3 * m_nGPNum, "wb");

	if (FIBER_INTER || FIBER_CLOSE) {
		count = 0;
		strcpy(file, path);
		strcat(file, "gp.fib");
		for (int j = 0; j < m_nGPNum; j++) {
			gp[count++] = m_spGPFib[j].x;
			gp[count++] = m_spGPFib[j].y;
			gp[count++] = m_spGPFib[j].z;
		}
		writeFile(file, gp, sizeof(double), 3 * m_nGPNum, "wb");
	}
	delete[] gp;

	strcpy(file, path);
	strcat(file, "gp.wei");
	writeFile(file, m_dpGPWei, sizeof(double), m_nGPNum, "wb");
	strcpy(file, path);
	strcat(file, "gp_inf.num");
	writeFile(file, m_npGPNBNum, sizeof(int), m_nGPNum, "wb");
	strcpy(file, path);
	strcat(file, "gp_inf.seq");
	writeFile(file, m_npGPNB, sizeof(int), m_nCurrentGP, "wb");
	strcpy(file, path);
	strcat(file, "gp.shape");
//	test("sha.bin",m_dpGPsha,sizeof(double),m_nCurrentGP,"wb");
	for (int j = 0; j < m_nGPNum; j++)
		m_mpGPsha[j].save(file, "a+b");
}

//*****************************************************************
//read in gpt info, when an existent set of gpts is used
//only cor, shaD, inf_num & inf_seq are needed for assembing while matrix
//*****************************************************************
void MFree::readGPInfo(char* path) {
	char* file = new char[MAX_FILENAME];
	int count = 0;

	double* gp = NULL;
	strcpy(file, path);
	strcat(file, "gp.cor");
	readFile(file, &gp, sizeof(double), 3, m_nGPNum);
	m_spGPCor = new cor[m_nGPNum];
	for (int j = 0; j < m_nGPNum; j++) {
		m_spGPCor[j].x = gp[count++];
		m_spGPCor[j].y = gp[count++];
		m_spGPCor[j].z = gp[count++];
	}
	delete[] gp;
	gp = NULL;

	strcpy(file, path);
	strcat(file, "gp_inf.num");
	readFile(file, &m_npGPNBNum, sizeof(int), 1, m_nGPNum);

	strcpy(file, path);
	strcat(file, "gp_inf.seq");
	readFile(file, &m_npGPNB, sizeof(int), 1, m_nCurrentGP);

	strcpy(file, path);
	strcat(file, "shaD.bin");
	readFile(file, &gp, sizeof(double), 3, m_nCurrentGP);
	m_mpGPshaD = new Matrix[m_nGPNum];
	count = 0;
	for (int j = 0; j < m_nGPNum; j++) {
		m_mpGPshaD[j].setData(&gp[count], 3, m_npGPNBNum[j]);
		count += 3 * m_npGPNBNum[j];
	}
	/*
	 count = 0;
	 strcpy(file,path);
	 strcat(file,"gp.fib");
	 readFil(file,&gp,sizeof(double),3,m_nGPNum);
	 for (int j=0; j<m_nGPNum; j++)
	 {
	 m_spGPFib[j].x = gp[count++];
	 m_spGPFib[j].y = gp[count++];
	 m_spGPFib[j].z    gp[count++];
	 }
	 delete []gp;

	 strcpy(file,path);
	 strcat(file,"gp.wei");
	 readFil(file,&gp,sizeof(double),3,m_nGPNum);
	 test(file,m_dpGPWei,sizeof(double),m_nGPNum,"wb");

	 */

}

//*******************************************************
//assTrans
//assemble transfer function got from each gp into global transfer function
//important function for the whole process of meshfree, integrate all the sub-functions together
//*******************************************************
void MFree::assTrans(BEM& bem, char* output) {
//	mwIndex col;
//	double* con = NULL;
	cor gpCor, gpWei;
	//first, genearate gp, wei, and corresponding shape functions
	//loop through all meshes
	clock_t s, f;
	double d;
	int j, k, l, m, n, o;

	cout << "Regenerate gauss points ? (y/n)" << endl;
	char c;
	cin >> c;
	if (c == 'y' || c == 'Y') {
		double* data = NULL;
		m_spGPCor = new cor[MAX_GP];
		m_dpGPWei = new double[MAX_GP];
		m_npGPNB = new int[MAX_GP * MAX_INF];  //sequnce of NB nodes for each gp
		m_npGPNBNum = new int[MAX_GP];  //number of NB nodes for each gp
		m_mpGPsha = new Matrix[MAX_GP];
		m_mpGPshaD = new Matrix[MAX_GP]; //[MAX_GP*MAX_INF*3];

		if (!FIBER_INTER && !FIBER_CLOSE) {
			int tempt;
			int pos = 0;
			readFile("../BEMForward/GP/gp_true.fib", &data, sizeof(double), 3,
					tempt);
			m_spGPFib = new cor[tempt];
			for (int i = 0; i < tempt; i++) {
				m_spGPFib[i].x = data[pos++];
				m_spGPFib[i].y = data[pos++];
				m_spGPFib[i].z = data[pos++];
			}
			if (data)
				delete[] data;
		} else
			m_spGPFib = new cor[MAX_GP];
		s = clock();
		cout
				<< "**********shape function construction*********************************"
				<< endl;

		m_ptrGP = new gp(MAX_GP);

		for (j = 0; j < m_npMesNum[1]; j++) {
			cout << "**********j =**************" << j << endl;
			for (k = 0; k < m_npMesNum[2]; k++) {
				cout << "**k =***" << k << endl;
				for (l = 0; l < m_npMesNum[3]; l++) {
					cout << l << endl;
					//info related to specific mesh
					m_ptrMes->setMesSeq(j, k, l, *this);
					//judge the status of the mesh
					//if totally outside, discard this mesh
					if (m_ptrMes->m_nSta != 0) {
						//loop through all the gps inside the mesh
						for (m = 0; m < m_ptrMes->m_nGPNum; m++) {
							for (n = 0; n < m_ptrMes->m_nGPNum; n++) {
								for (o = 0; o < m_ptrMes->m_nGPNum; o++) {
									//get info about specific gp
									// cout<<"generating global gauss pts"<<endl;
									m_ptrMes->genGloGP(gpCor, gpWei, m, n, o);
									//update gp object
									if (!isout(gpCor)) {
										//cout<<"re-initializing global gauss pts"<<endl;
										m_ptrGP->reIni(gpCor, gpWei, *this);

										//cout<<"calculating single transfer matrix"<<endl;
										m_ptrGP->sinTrans(*this);
									} //each in gpt
								}    //z
							}  //y
						}   //x
					}  //each bg box partially/fully inside
				}  //z
			}  //y
		}  //x
		//	write out gp nodes for visualization***********************************
		cout << "number of gs points used: " << m_nGPNum << endl;
		//writeGPInfo("/Users/maomaowlw/Research/Mywork/Workspace/BEMForward/GP/");
		delete m_ptrGP;
	} else {
		readGPInfo("../BEMForward/GP/");
		cout << "number of gpts used: " << m_nGPNum << endl;
	}

	f = clock();
	d = double(f - s) / CLOCKS_PER_SEC;
	::cout << "time for shape functions" << d << endl;
	//begin to assemble the whole matrix******************************************
	s = clock();
	double r, dx, dy, dz;
	int pos_sur, pos_bem, pos_gp;
	double* ptr = NULL;
	cout
			<< "*****transfer matrix assemble**************************************"
			<< endl;
	for (j = 1; j <= bem.m_nNodNum; j++) {
		if (!bem.isForbidden(j)) {
			cerr << j;
			//info related to specific surface point
			setSurPoiInf(j, d2c(&bem.m_dpNodCor[3 * (j - 1)])); //,d2c(bem.m_dpnodCon[j-1]));
			//start from the first gp
			pos_sur = 0;
			pos_bem = (j - 1) * m_nNodNum;
			for (l = 0; l < m_nGPNum; l++) {
				//construct [rx,ry,rz]/(4*pi*r^2)
				dx = m_spGPCor[l].x - m_sSurPoiCor.x;
				dy = m_spGPCor[l].y - m_sSurPoiCor.y;
				dz = m_spGPCor[l].z - m_sSurPoiCor.z;
				r = sqrt(dx * dx + dy * dy + dz * dz);
				r = 1 / (4 * PI * r * r * r * Deff);
				m_dpR[0] = dx * r;
				m_dpR[1] = dy * r;
				m_dpR[2] = dz * r;
				//read in shape functions
				ptr = m_mpGPshaD[l].getDataPtr();
				pos_gp = 0;
				//assemble into global matrix
				for (k = 0; k < m_npGPNBNum[l]; k++) {
					m_dpTrans[pos_bem + m_npGPNB[pos_sur] - 1] += m_dpR[0]
							* ptr[pos_gp] + m_dpR[1] * ptr[pos_gp + 1]
							+ m_dpR[2] * ptr[pos_gp + 2];
					pos_sur++;
					pos_gp += 3;
				}
			}
		}
	}
	pos_bem = m_nNodNum * bem.m_nNodNum;
	for (j = 0; j < m_nNodNum; j++)
		m_dpTrans[pos_bem++] = 0;

//	test("MFree_original.bin",m_dpTrans,sizeof(double),m_nNodNum*(bem.m_nNodNum+1),"wb");

	rearrange(bem);

	writeFile(output, m_dpTrans_arranged, sizeof(double),
			m_nNodNum * (bem.m_nNodNum + 1 - bem.m_nForNum), "wb");
	f = clock();
	d = double(f - s) / CLOCKS_PER_SEC;
	::cout << "time for assemble" << d << endl;

	//Div(u) on x, y, z directions, on gp points
	/*	double* dx_gp = new double[m_nGPNum*m_nNodNum];
	 double* dy_gp = new double[m_nGPNum*m_nNodNum];
	 double* dz_gp = new double[m_nGPNum*m_nNodNum];
	 double* ptr;
	 int seq;
	 pos=0;
	 for( l=0; l<m_nGPNum; l++)
	 {
	 ptr = mxGetPr(m_mpGPshaD[l].GetData());
	 for(j=0; j<m_npGPNBNum[l]; j++)
	 {
	 seq = l*m_nNodNum+m_npGPNB[pos]-1;
	 dx_gp[seq] = *ptr++;
	 dy_gp[seq] = *ptr++;
	 dz_gp[seq] = *ptr++;
	 pos++;
	 }
	 }
	 wriFil(".\\heart836torso700\\test\\Dux.bin",dx_gp,sizeof(double),m_nGPNum*m_nNodNum,"wb");
	 wriFil(".\\heart836torso700\\test\\Duy.bin",dy_gp,sizeof(double),m_nGPNum*m_nNodNum,"wb");
	 wriFil(".\\heart836torso700\\test\\Duz.bin",dz_gp,sizeof(double),m_nGPNum*m_nNodNum,"wb");

	 delete []dx_gp;
	 delete []dy_gp;
	 delete []dz_gp;
	 */
}

//**********************************************************
//rearrange: get rid of elements corresponding to inner nodes in BEM
//**********************************************************
void MFree::rearrange(BEM& bem) {
	if (FORBIDDEN) {
		int gpos, pos = 0;
		for (int i = 0; i < bem.m_nNodNum; i++) {
			if (!bem.isForbidden(i + 1)) {
				gpos = i * m_nNodNum;
				for (int j = 0; j < m_nNodNum; j++)
					m_dpTrans_arranged[pos++] = m_dpTrans[gpos + j];
			}
		}
		for (int k = 0; k < m_nNodNum; k++)
			m_dpTrans_arranged[pos++] = m_dpTrans[gpos++];

		if (pos != m_nNodNum * (bem.m_nNodNum + 1 - bem.m_nForNum)) {
			cout << "errors during rearrange MFree" << endl;
			exit(0);
		}
	} else {
		for (int i = 0; i < (bem.m_nNodNum + 1) * m_nNodNum; i++)
			m_dpTrans_arranged[i] = m_dpTrans[i];
	}
}

