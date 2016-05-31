/*
 * gp.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: jxx2144
 */

#include "inclu.h"
//#include"BEMCell.h"
#include "external.h"
#include "gp.h"
//#include "BEM.h"
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
bgMesh::bgMesh(MFree& mfree)
{
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
bgMesh::~bgMesh()
{

	delete []m_dpNatGPCor;
	delete []m_dpNatGPWei;
}

//*******************************************************
//setMesSeq
//change the mesh domain as loop continues, i.e. to sweep through all the meshes
//*******************************************************
void bgMesh::setMesSeq(int i,int j,int k, MFree& mfree)
{
		// generate domain of specified mesh
	xl = mfree.m_spBouCor[0].x + i*mfree.m_nMesInt;
	xu = xl + mfree.m_nMesInt;

	yl = mfree.m_spBouCor[0].y + j*mfree.m_nMesInt;
	yu = yl + mfree.m_nMesInt;

	zl = mfree.m_spBouCor[0].z + k*mfree.m_nMesInt;
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
void bgMesh::comSta(MFree& mfree)
{
	m_nSta = 2;
	m_nGPNum = 0; // re-initialize

	for(int i=0; i<mfree.m_nNodNum; i++)
	{
		//at least one nodal is in the mesh, so it's not totally out
		if( mfree.m_spNodCor[i].x >= xl && mfree.m_spNodCor[i].x <= xu
		&& mfree.m_spNodCor[i].y >= yl && mfree.m_spNodCor[i].y <= yu
		&& mfree.m_spNodCor[i].z >= zl && mfree.m_spNodCor[i].z <= zu)

			m_nGPNum = m_nGPNum + 1;
	}

	if(m_nGPNum==0)
	  m_nSta = 0;

		// generate natural gp coordinates&weights, just for 1D info

	else
	{
		m_nGPNum = floor (pow(m_nGPNum, 1.0/3.0 ) + 2);
		for(int j=0; j<MAX_NAT_GP; j++)
		{
			m_dpNatGPCor[j] = 0;
			m_dpNatGPWei[j] = 0;
		}
		switch(m_nGPNum)
		{
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
void bgMesh::genGloGP(cor &gpcor,cor& gpwei, int i, int j,int k)
{
	// x = ( xl + xu) /2 + gp*(xu-xl)/2
	gpcor.x = ( xl + xu )/2.0 + m_dpNatGPCor[i]*( xu - xl)/2.0;
	gpcor.y = ( yl + yu )/2.0 + m_dpNatGPCor[j]*( yu - yl)/2.0;
	gpcor.z = ( zl + zu )/2.0 + m_dpNatGPCor[k]*( zu - zl)/2.0;
	// neww = (xu-xl)*oldw/2
	gpwei.x = m_dpNatGPWei[i]*(xu-xl)/2.0;
	gpwei.y = m_dpNatGPWei[j]*(yu-yl)/2.0;
	gpwei.z = m_dpNatGPWei[k]*(zu-zl)/2.0;
}


//**************************************************************************************************************-
//class:nodal
//**************************************************************************************************************
//*******************************************************
//ctors
//get the cor of the nodal according to input seq, and generate p for computing shape function
//*******************************************************
nodal::nodal():p(4,1)
{
//	p = new double[4];
         // no assignment, just generate the structure
}

//*******************************************************
//dtor
//*******************************************************
nodal::~nodal()
{
//	delete p;
}

//*******************************************************
//reIni
//get nodal info for each nodal according to its sequence
//*******************************************************-------
void nodal::reIni(int seq, MFree& mfree)
{
	m_cor = mfree.m_spNodCor[seq-1];

//	p = vertcat(1, m_cor.x, m_cor.y, m_cor.z,m_cor.x*m_cor.y, m_cor.x*m_cor.z, m_cor.y*m_cor.z,
//	m_cor.x*m_cor.x, m_cor.y*m_cor.y, m_cor.z*m_cor.z);
	p.setData(1,0,0);
	p.setData(m_cor.x,1,0);
	p.setData(m_cor.y,2,0);
	p.setData(m_cor.z,3,0);
}


//****************************************************************************************************************
//class applicable for Poly (p,p_dx,p_dy,p_dz); Shape; A, B
//***************************************************************************************************************
//without dimension information, intialize as empty matrix
MatrixGroup::MatrixGroup():d(),dx(),dy(),dz()
{}

MatrixGroup::MatrixGroup(int row, int col):d(row,col),dx(row,col),dy(row,col),dz(row,col)
{}

void MatrixGroup::Initialize(int row, int col)
{
	d.setData(0.0,row,col);
	dx.setData(0.0,row,col);
	dy.setData(0.0,row,col);
	dz.setData(0.0,row,col);
}

void MatrixGroup::reSet()
{
	d.reSet(0.0);
	dx.reSet(0.0);
	dy.reSet(0.0);
	dz.reSet(0.0);
}


MatrixGroup::~MatrixGroup()
{}
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
gp::gp():D0(3,3),D(3,3),Q(3,3),shape(),p(4,1),sha(),A(4,4),B(),invA(4,4)
{
	Q.setData(1.0,0,0);
	Q.setData(1.0,1,1);
	Q.setData(1.0,2,2);
	if(ISO == false)
	{
		D0.setData(DI_S,0,0);
		D0.setData(DI_S,1,1);
		D0.setData(DI_L,2,2);
	}
	else
	{
		D0.setData(DI,0,0);
		D0.setData(DI,1,1);
		D0.setData(DI,2,2);
	}

	m_ptrNod = new nodal();

	m_nSur = 0;
	m_npSurSeq = NULL;
	m_nodWei = NULL;
	m_nodWeiDx= NULL;
	m_nodWeiDy = NULL;
	m_nodWeiDz = NULL;

	m_nCurrent = 0;
}

//****************************************************
//for a bunch of gpts, use fixed length to initialize all memory related to
//influence nodes (use MAX_INF)
//****************************************************
gp::gp(int num):D0(3,3),D(3,3),Q(3,3),shape(3,num),p(4,1),sha(1, num),A(4,4),B(4,num),invA(4,4)
{
	Q.setData(1.0,0,0);
	Q.setData(1.0,1,1);
	Q.setData(1.0,2,2);
	if(ISO == false)
	{
		D0.setData(DI_S,0,0);
		D0.setData(DI_S,1,1);
		D0.setData(DI_L,2,2);
	}
	else
	{
		D0.setData(DI,0,0);
		D0.setData(DI,1,1);
		D0.setData(DI,2,2);
	}

	m_ptrNod = new nodal();

	m_nSur = 0;
	m_npSurSeq = new int[num];
	m_nodWei = new double[num];
	m_nodWeiDx= new double[num];
	m_nodWeiDy = new double[num];
	m_nodWeiDz = new double[num];

	m_nCurrent = 0;
}

//*******************************************************
//dtor
//*******************************************************
gp::~gp()
{
	delete []m_npSurSeq;
	delete []m_nodWei;
	delete []m_nodWeiDx;
	delete []m_nodWeiDy;
	delete []m_nodWeiDz;

	//delete m_ptrNod;
}

//*******************************************************
//reIni
//re-initialize, i.e. re-set the info of gp when loop through each gp
//*******************************************************
void gp::reIni(cor& gpCor, cor& gpWei, MFree& mfree)
{
	m_gpCor = gpCor;
	m_dWei = gpWei.x*gpWei.y*gpWei.z;

	m_nCurrent++;

	Poly3D(p.d,m_gpCor.x, m_gpCor.y, m_gpCor.z);
	Poly3D_dx(p.dx,m_gpCor.x, m_gpCor.y, m_gpCor.z);
	Poly3D_dy(p.dy,m_gpCor.x, m_gpCor.y, m_gpCor.z);
	Poly3D_dz(p.dz,m_gpCor.x, m_gpCor.y, m_gpCor.z);

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
void gp::getSur(MFree& mfree)
{
	double a, r, dx, dy, dz;
	m_nSur = 0;
	for(int k=0; k<MAX_INF; k++)
	{
		m_npSurSeq[k] = 0;
		m_nodWei[k] = 0;
		m_nodWeiDx[k] = 0;
		m_nodWeiDy[k] = 0;
		m_nodWeiDz[k] = 0;
	}
	//loop through each nodal, to check which ones influence domain cover gp
	for ( int i = 0; i< mfree.m_nNodNum; i++)
	{  //the distance between gp and current nodals
		a = DMAX*mfree.m_dpNodDom[i];
		dx = mfree.m_spNodCor[i].x - m_gpCor.x;
		dy = mfree.m_spNodCor[i].y - m_gpCor.y;
		dz = mfree.m_spNodCor[i].z - m_gpCor.z;
		r = sqrt(dx*dx + dy*dy + dz*dz);

		if (r > 0.00001)
		{
			dx /= (a*r);
			dy /= (a*r);
			dz /= (a*r);
			r /= a;

			if( r <= 0.5 )
    		{
    			if(m_nSur+1>MAX_INF)
    			{
    				cout<<"influencing nodes out of limits"<<endl;
    		        exit(0);
    			}
    			else
    			{
			// w = 2/3 - 4r^2 + 4r^3
    				m_nodWei[m_nSur] = 2.0/3.0 - 4*r*r + 4*r*r*r;
		    // wx = -8rx + 12r*rx    ( rx = (x - xi)/(a^2*r))
    				m_nodWeiDx[m_nSur] = (-8*r + 12*r*r)*dx;
    		        m_nodWeiDy[m_nSur] = (-8*r + 12*r*r)*dy;
    		        m_nodWeiDz[m_nSur] = (-8*r + 12*r*r)*dz;
		// wxx = -8 + 12(rx^2/r + r)
//		m_nodWeiDD.x = (-8+24*r)*dx*dx + m_nodWeiD.x*ddx/dx;
//		m_nodWeiDD.y = (-8+24*r)*dy*dy + m_nodWeiD.x*ddy/dy;
//		m_nodWeiDD.z = (-8+24*r)*dz*dz + m_nodWeiD.z*ddz/dz;
    		        m_npSurSeq[m_nSur++] = i+1;
    		     }
    	      }
         	  else if( 0.5<r && r<=1 )
    		{
    			if(m_nSur+1>MAX_INF)
    			{
    				cout<<"influencing nodes out of limits"<<endl;
    				wait();
    		        exit(0);
    			}
    			else
    			{
			   // w = 4/3 - 4r + 4r^2 - 4/3r^3
	    			m_nodWei[m_nSur] = 4.0/3.0 - 4*r + 4*r*r - 4.0/3.0*r*r*r;
		      // wx = -4rx/r + 8rx - 4rrx
                    m_nodWeiDx[m_nSur] = (-4+8*r-4*r*r)*dx;
    		        m_nodWeiDy[m_nSur] = (-4+8*r-4*r*r)*dy;
    		        m_nodWeiDz[m_nSur] = (-4+8*r-4*r*r)*dz;
		      //wxx = 8 -4(1/r - rx^2/r^3) - 4(rx^2 + r)
//		m_nodWeiDD.x = (8 -8*r)*dx*dx + m_nodWeiD.x*ddx/dx;
//		m_nodWeiDD.y = (8 -8*r)*dy*dy + m_nodWeiD.y*ddy/dy;
//		m_nodWeiDD.z = (8 -8*r)*dz*dz + m_nodWeiD.z*ddz/dz;
    			    m_npSurSeq[m_nSur++] = i+1;
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
void gp::sinTrans(MFree& mfree)
{
	Matrix tmp(4,4);

	//loop through every supportive point
	//construct A,B, and their derivatives
	for ( int i = 0; i< m_nSur; i++)
	{
		//compute A,B and their derivatives nodal by nodal
		m_ptrNod->reIni(m_npSurSeq[i], mfree);
		tmp = m_ptrNod->p*(~(m_ptrNod->p));
		A.d += tmp*m_nodWei[i];
		A.dx += tmp*m_nodWeiDx[i];
		A.dy += tmp*m_nodWeiDy[i];
		A.dz += tmp*m_nodWeiDz[i];

        B.d.setSub((m_ptrNod->p)*m_nodWei[i],0,3,i,i);
		B.dx.setSub((m_ptrNod->p)*m_nodWeiDx[i],0,3,i,i);
		B.dy.setSub((m_ptrNod->p)*m_nodWeiDy[i],0,3,i,i);
		B.dz.setSub((m_ptrNod->p)*m_nodWeiDz[i],0,3,i,i);
	}
	//compute ddSha
//	cout<<"A inverse"<<endl;

	invA.d = A.d.inv();

//	cout<<"sha.d"<<endl;
	sha.d.setSub((~(p.d))*invA.d*B.d.sub(0,3,0,m_nSur-1), 0, 0, 0, m_nSur-1);
//	sha.d.save("../BEMForward/GP/sha.bin","a+b");
	if (checkNan(sha.d))
	{
		cout<<"Nan happens @ sha"<<endl;
		wait();
	}

	int dsum = sha.d.sum() - 1.0;
	if( dsum > 0.00001|| dsum < -0.00001)
	{
		cout<<"Error: sum of shape is not equal to 1 on gpt"<<endl;
		wait();
	}

	invA.dx = - invA.d*A.dx*invA.d;
//	invAdd[0] = -(invAd[0]*Ad[0]*invA + invA*Add[0]*invA + invA*Ad[0]*invAd[0]);
	invA.dy = - invA.d*A.dy*invA.d;
//	invAdd[1] = -(invAd[1]*Ad[1]*invA + invA*Add[1]*invA + invA*Ad[1]*invAd[1]);
	invA.dz = - invA.d*A.dz*invA.d;
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
//	cout<<"sha.dxyz"<<endl;
	sha.dx.setSub((~p.dx)*invA.d*B.d.sub(0,3,0,m_nSur-1) + (~p.d)*invA.dx*B.d.sub(0,3,0,m_nSur-1) + (~p.d)*invA.d*B.dx.sub(0,3,0,m_nSur-1),0,0,0,m_nSur-1);
	sha.dy.setSub((~p.dy)*invA.d*B.d.sub(0,3,0,m_nSur-1) + (~p.d)*invA.dy*B.d.sub(0,3,0,m_nSur-1) + (~p.d)*invA.d*B.dy.sub(0,3,0,m_nSur-1),0,0,0,m_nSur-1);
	sha.dz.setSub((~p.dz)*invA.d*B.d.sub(0,3,0,m_nSur-1)+ (~p.d)*invA.dz*B.d.sub(0,3,0,m_nSur-1) + (~p.d)*invA.d*B.dz.sub(0,3,0,m_nSur-1),0,0,0,m_nSur-1);

	if (checkNan(sha.dx))
	{
		cout<<"Nan happens @ sha_dx"<<endl;
		wait();
	}
	if (checkNan(sha.dy))
	{
		cout<<"Nan happens @ sha_dy"<<endl;
		wait();
	}
	if (checkNan(sha.dz))
	{
		cout<<"Nan happens @ sha_dz"<<endl;
		wait();
	}

//	cout<<"set shape function"<<endl;
	//if(ISO==0)
	//	setGPCon(mfree,NULL);
	shape.setSub(sha.dx.sub(0,0,0,m_nSur-1), 0,0,0,m_nSur-1);
	shape.setSub(sha.dy.sub(0,0,0,m_nSur-1),1,1,0,m_nSur-1);
	shape.setSub(sha.dz.sub(0,0,0,m_nSur-1),2,2,0,m_nSur-1);

//	mwArray shape = vertcat(ddSha[0], ddSha[1], ddSha[2]);
//	double* ptr = shape.getDataPtr();
//	writeFile("../BEMForward/GP/shapeD.bin",ptr,sizeof(double),3*m_nSur,"a+b");

   //***store into mfree class************
//	cout<<"store gp"<<endl;

	if(mfree.m_nGPNum+1>MAX_GP)
	{
		cout<<"not enough memeory for gps: "<<endl;
		exit(0);
	}
	mfree.m_spGPCor[mfree.m_nGPNum] = m_gpCor; // gp cor,record after first sweep
	//mfree.m_spGPFib[mfree.m_nGPNum] = m_gpFib;
	mfree.m_dpGPWei[mfree.m_nGPNum] = m_dWei;
	mfree.m_npGPNBNum[mfree.m_nGPNum] = m_nSur;  //number of NB nodes for each gp
	mfree.m_mpGPshaD[mfree.m_nGPNum].setData(0.0,3,mfree.m_nNodNum);
	mfree.m_mpGPshaD[mfree.m_nGPNum].setSub((shape.sub(0,2,0,m_nSur-1)),0,2,0,m_nSur-1);
	//mfree.m_mpGPshaD->setSubID((shape.sub(0,2,0,m_nSur-1)),3*mfree.m_nGPNum,3*mfree.m_nGPNum+2,m_nSur, m_npSurSeq);

	int gg = mfree.m_nGPNum * MAX_INF;
	mfree.m_nGPNum++;
	for(int j=0; j<m_nSur; j++)
	{
	    mfree.m_npGPNB[gg] = m_npSurSeq[j];      //sequnce of NB nodes for each gp
//		mfree.m_dpGPsha[mfree.m_nCurrentGP] = sha.d(0,j);
		gg++;
	}

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
MFree::MFree(char* file)
{
	//read in data-------------------------------------------------------------
	//nodal cor

	//open file, and get the size of the file, then correpondingly # of data and further # of nodals
	//nodal cor
    double* data = NULL;
	int i,pos=0;
	//char* file = new char[MAX_FILENAME];
	readFile(file,&data,sizeof(double),3,m_nNodNum);
	m_spNodCor = new cor[m_nNodNum];
	//read in nodal by nodal
	for(i=0; i<m_nNodNum; i++)
	{
	    m_spNodCor[i].x = data[pos++];
	    m_spNodCor[i].y = data[pos++];
	    m_spNodCor[i].z = data[pos++];
	}
	if(data)
		delete []data;
	data = NULL;


	//*****************************
	//***GP info allocated later, dependent on whether generate new or not****************
    m_spGPCor = NULL;
	m_dpGPWei = NULL;
	m_nGPNum = 0;
	m_nCurrentGP = 0;
	m_npGPNB = NULL;      //sequnce of NB nodes for each gp
	m_npGPNBNum = NULL;  //number of NB nodes for each gp
	m_mpGPshaD = NULL; //[MAX_GP*MAX_INF*3];
	m_ptrGP = NULL;
//	m_dpGPsha = new double[MAX_GP*MAX_INF];

	m_dpR = new double[3];


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
MFree::~MFree()
{
	delete []m_spNodCor;

	delete []m_dpNodDom;

    delete []m_spGPCor;

    if (m_dpGPWei)
		delete []m_dpGPWei;
	delete []m_npGPNB;
	delete []m_npGPNBNum;
//	delete []m_dpGPsha;
	delete []m_mpGPshaD;

	delete []m_spBouCor;
	delete []m_npMesNum;
	delete []m_dpR;

}


//*******************************************************
//getNodInfDom:
//according to specified # of neighbours for each nodal
//*******************************************************
void MFree::getNodInfDom(int num)
{
	m_dpNodDom = new double[m_nNodNum];
	double* inf = new double[num];       // record of influence neighbors in ascent sequence
	double r=0;
	for ( int l=0; l<m_nNodNum; l++)
	{
		for(int i=0; i<num; i++)
			inf[i] = MAX_R;
		//loop through every nodal
		for(int j=0; j<m_nNodNum; j++)
		{
			//loop through every nodal
			if( l!=j)
			{
				r = sqrt(power(m_spNodCor[l].x - m_spNodCor[j].x) + power(m_spNodCor[l].y - m_spNodCor[j].y)
				    + power(m_spNodCor[l].z - m_spNodCor[j].z));
				//check through the record queue, and add/update accordingly
				//determine the pos in record queue
				int k = 0;
				for( k=0; k<num && r>inf[k]; k++)
				{}
				//update record queue( insert in new neighbor at correct pos, and move following neighbors backward)
				if (k < num)
				{
					for( int m= (num-1); m>k; m--)
						inf[m] = inf[m-1];
					inf[k] = r;
				}
			}
		}
		//get radius for influence domain
		//m_dpNodDom[l] = inf[num-1];
		m_dpNodDom[l] = 0;
		for (int i = 0; i<num; i++)
			m_dpNodDom[l] += inf[i];
		m_dpNodDom[l] = m_dpNodDom[l]/num;

	}
	delete []inf;
}



//*******************************************************
//setSurPoiInf
//set the current surface point info, which is got from external application, here BEM object
//*******************************************************
void MFree::setSurPoiInf(int seq, cor poiCor)//, cor poiCon)
{
	m_nSurPoiSeq = seq;
	m_sSurPoiCor = poiCor;
//	m_dSurPoiCon = poiCon.x;
}

//*******************************************************
//genBgMes
//generate bg mesh according to read-in data of nodals, actually target is to get the boundary cor for mesh, and
//basic info about mesh, eg. #gp, #meshes
//*******************************************************
void MFree::genBgMes()
{
	//pre-allocate storage
	m_spBouCor = new cor[2];

	m_spBouCor[0] = m_spNodCor[0];
	m_spBouCor[1] = m_spNodCor[0];

	//loop through every tetra nodes to determine the lower bound and upper bound of bg mesh
	for(int i=0; i< m_nTetNodNum; i++)
	{
		if(m_spBouCor[0].x > m_spTetCor[i].x)
			m_spBouCor[0].x = m_spTetCor[i].x;
		if(m_spBouCor[0].y > m_spTetCor[i].y)
			m_spBouCor[0].y = m_spTetCor[i].y;
		if(m_spBouCor[0].z > m_spTetCor[i].z)
			m_spBouCor[0].z = m_spTetCor[i].z;
		if(m_spBouCor[1].x < m_spTetCor[i].x)
			m_spBouCor[1].x = m_spTetCor[i].x;
		if(m_spBouCor[1].y < m_spTetCor[i].y)
			m_spBouCor[1].y = m_spTetCor[i].y;
		if(m_spBouCor[1].z < m_spTetCor[i].z)
			m_spBouCor[1].z = m_spTetCor[i].z;
	}
	//compute the total# mesh
	//get the directional #  -- round off into higher integral
	double dx = (m_spBouCor[1].x - m_spBouCor[0].x);
	double dy = (m_spBouCor[1].y - m_spBouCor[0].y);
	double dz = (m_spBouCor[1].z - m_spBouCor[0].z);
	m_npMesNum[1] = ceil(dx/ m_nMesInt);
	m_npMesNum[2] = ceil(dy/ m_nMesInt);
	m_npMesNum[3] = ceil(dz/ m_nMesInt);
	// adjust the boundary cor to make integral# mesh
	dx = m_npMesNum[1]*m_nMesInt - dx;
	dy = m_npMesNum[2]*m_nMesInt - dy;
	dz = m_npMesNum[3]*m_nMesInt - dz;
	m_spBouCor[0].x -= dx/2;
	m_spBouCor[0].y -= dy/2;
	m_spBouCor[0].z -= dz/2;
	m_spBouCor[1].x += dx/2;
	m_spBouCor[1].y += dy/2;
	m_spBouCor[1].z += dz/2;
	//total #
	m_npMesNum[0] = m_npMesNum[1]*m_npMesNum[2]*m_npMesNum[3];
    cout<<"Number of background meshes: "<<m_npMesNum[0];
	//total # of gp
	//m_nGPTot = power(m_nGPNum)*m_nGPNum * m_npMesNum[0];
}
//*******************************************************
//assTrans
//assemble transfer function got from each gp into global transfer function
//important function for the whole process of meshfree, integrate all the sub-functions together
//*******************************************************
void MFree::assTrans(char* file)
{
	cor gpCor, gpWei;
	//first, genearate gp, wei, and corresponding shape functions
	//loop through all meshes


	int j,lll, jj, k,l;

	//cout<<"Regenerate gauss points ? (y/n)"<<endl;
	char c;
	//cin>>c;
	/////////read guess point
	//char* file = new char[MAX_FILENAME];
	int count = 0;
	double* gpt = NULL;

	//strcat(file,"cube1.cor");
	readFile(file,&gpt,sizeof(double),3,lll);
	m_spGPCor = new cor[lll];
	for (int j=0; j<lll; j++)
	{
		m_spGPCor[j].x = gpt[count++];
		m_spGPCor[j].y = gpt[count++];
		m_spGPCor[j].z = gpt[count++];
	}
    delete []gpt;
	////////////////////////
	c = 'y';
	if (c == 'y' || c == 'Y')
	{
		double* data = NULL;
		//m_spGPCor= new cor[MAX_GP];
    	m_dpGPWei = new double[lll];
    	m_npGPNB = new int[lll*MAX_INF];      //sequnce of NB nodes for each gp
    	m_npGPNBNum = new int[lll];  //number of NB nodes for each gp
    	m_mpGPshaD = new Matrix[lll]; //[MAX_GP*MAX_INF*3];


		//s = clock();
    	cout<<"**********shape function construction*********************************"<<endl;

		m_ptrGP = new gp(lll);
	/*	for( j=0; j<m_npMesNum[1]; j++)
				{
		        	cout<<"**********j =**************"<<j<<endl;
		    		for( k=0; k<m_npMesNum[2]; k++)
					{
		    			cout<<"**k =***"<<k<<endl;
		    		    for( l=0; l<m_npMesNum[3]; l++)
						{
		    		    	cout<<l<<endl;
		    			    //info related to specific mesh
		    	    	    m_ptrMes->setMesSeq(j,k,l,*this);
		    			    //judge the status of the mesh
		    			    //if totally outside, discard this mesh
		     				if(m_ptrMes->m_nSta != 0)
							{

			    				    	  //get info about specific gp
										  //cout<<"global gp"<<endl;
		     								int m = 0;
		     								int n = 0;
		     								int o = 0;
			    					      m_ptrMes->genGloGP(gpCor,gpWei,m,n,o);
*/
											for( jj=0; jj<lll; jj++)
											{
												//cout<<"**********j =**************"<<j<<endl;
													gpCor.x = m_spGPCor[jj].x;
													gpCor.y = m_spGPCor[jj].y;
													gpCor.z = m_spGPCor[jj].z;

													 m_ptrGP->reIni(gpCor,gpWei,*this);
													 m_ptrGP->sinTrans(*this);
													// m_mpGPshaD->setSubID((m_ptrGP->shape.sub(0,2,0,m_npGPNBNum[j]-1)),3*j,3*j+2,m_npGPNBNum[j], m_npGPNB);
													 ::cout<<jj<<"shape: "<<m_ptrGP->shape.sub(0,2,0,m_npGPNBNum[jj]-1).getDataPtr()[15]<<endl;


											}  //x
					//		}// if
				//		}// l
				//	}//k
			//	}//j
       //	write out gp nodes for visualization***********************************
    	::cout<<"number of gs points used: "<<m_nGPNum<<endl;
    	delete m_ptrGP;
	}

//	f = clock();
	//d = double(f-s)/CLOCKS_PER_SEC;
	::cout<<"time for shape functions" ;//<< d<<endl;
    //begin to assemble the whole matrix******************************************
	//s = clock();


}



