// MFree.h: declaration of classes needed in MFree method

#ifndef _MFREE__H_
#define __MFREE__H_

//************************************************************************************************
//inclusions and references
//************************************************************************************************
#include "Matrix.h"

class MFree;

//************************************************************************************************
//declaratin of classes
//************************************************************************************************
//*******************************************************
//class: bgMesh
// used to specify each bgmesh according to its global sequence(coordinate of boundary)
// generate according gp; judge whether outof original volume domain
//*******************************************************
class bgMesh {
	friend class MFree;
	//functions
	//ctors/dtor
public:
	bgMesh(MFree& mfree);
	~bgMesh();

private:
	void setMesSeq(int, int, int, MFree& mfree);
	void genGloGP(cor& gpcor, cor& gpwei, int, int, int);
	void comSta(MFree& mfree);     // 0, totally out; 1, partially out; else, in

	//member
	//individual ID,specified each time
	double xl;
	double xu;
	double yl;
	double yu;
	double zl;
	double zu;
	int m_nSta;
	// shared info, only initialized once
	double* m_dpNatGPCor;
	double* m_dpNatGPWei;
	int m_nGPNum;         // # gp for each mesh
};

//*******************************************************
//class:nodal
//spedify each nodal: determine P matrix
//*******************************************************
class nodal {
	friend class gp;
	//ctors/dtor
public:
	nodal();
	~nodal();
	//re-intialization
private:
	void reIni(int, MFree&);

	//member
private:
	cor m_cor;
	Matrix p;   //    double* p;
};

//*******************************************************************************
//general class for shape, A, B which has **, **_dx, **_dy, **_dz
//********************************************************************************
class MatrixGroup {
public:
	MatrixGroup();
	MatrixGroup(int row, int col);
	~MatrixGroup();

	void Initialize(int row, int col);
	void reSet();
public:
	Matrix d;
	Matrix dx;
	Matrix dy;
	Matrix dz;
};

//*******************************************************
//class: gp
//used to specify each global gp(cor&wei); get it supportive nodals(seq and #); compute
// f(x) for this gp; assemble it into global transfer matrix
//*******************************************************
class gp {
	friend class nodal;
	friend class MFree;

	//ctors/dtor
public:
	gp();
	gp(int num);
	~gp();
private:
	//functions
	void reIni(cor&, double, MFree&);
	void reIni(cor&, cor&, MFree&);            // set the properties of each gp
	void setGPCon(MFree&, double*);
	void computeQ();
	void getSur(MFree&);
	void poly3D(Matrix* m, cor c);

	void sinTrans(MFree&);

	//members
	cor m_gpCor;
	cor m_gpFib;
	double m_dWei;
	//information for support nodes
	double* m_nodWei;
	double* m_nodWeiDx;
	double* m_nodWeiDy;
	double* m_nodWeiDz;
	int m_nSur;
	int* m_npSurSeq;         // supportive nodes and its sequences
	int m_nCurrent;

	Matrix D0;
	Matrix D;
	Matrix Q;
	Matrix shape;
	MatrixGroup p;
	MatrixGroup sha;
	MatrixGroup A;
	MatrixGroup B;
	MatrixGroup invA;

	nodal* m_ptrNod;      //pointer to a nodal
};

//*******************************************************
//class:MFree
// used to manage objects of gp,nodal & bgMesh to generate the whole transfer matrix
// read in data
// output: transfer matrx ; get input and compute body surface potential
//*******************************************************
class MFree {
	friend class gp;
	friend class bgMesh;
	friend class nodal;
	friend class BEM;
	friend class BEM_source;
	friend class BEM_nosource;

	//ctors/dtor
public:
	MFree(double, char**, BEM& bem); // specify gp in each mesh, mesh interval and input file
	~MFree();
private:
	//functions
	void genBgMes();                    //apply when initialize
	void getNodInfDom(int); //generate influence domain diameter for every nodal,apply when initialize
	void setSurPoiInf(int, cor); //set surface points info, apply for each surface points
	bool isout(cor gpCor); //detect whether gp( with cor gpCor) is in the geometry
	void rearrange(BEM& bem);

	void writeGPInfo(char* path);
	void readGPInfo(char* path);     //all gpt info saved under //GP// directory

public:
	void assTrans(BEM&, char* output); //assemble transfer matrix ( major function)

	//member
	//read-in data
private:
	cor* m_spNodCor;
	cor* m_spNodCon;      // nodal info, fiber orientation

	cor* m_spGPCor;
	cor* m_spGPFib;       //fib orientation of each gps.
	double* m_dpGPWei;     // gp cor and weis,record after first sweep
	int* m_npGPNB;      //sequnce of NB nodes for each gp
	int* m_npGPNBNum;  //number of NB nodes for each gp
//	double* m_dpGPsha;
	Matrix* m_mpGPsha;
	Matrix* m_mpGPshaD;
//	double* m_mpGPshaD;
	int m_nGPNum;
	int m_nCurrentGP;

	//variables computed in " getNodInfDom"
	double* m_dpNodDom;   //influence domain diameter for every nodal    
	// variables computed in" setSurPoiInf"
	int m_nSurPoiSeq; // sequence of surface point, used to determine which row of transfer matrix is computed
	cor m_sSurPoiCor; // coordinate of surface point, used in computation of f(x)
	double* m_dpR;
//	double m_dSurPoiCon;
	//variables computed in "genBgMes"
	double m_nMesInt;        // bg_mesh interval
	cor* m_spBouCor; // apex of boundary mesh: [(xl,yl,zl),(xu,yu,zu)]not correpsonds to real point coordinates
	int* m_npMesNum;     //#  mesh  ( total# , # in x direction, # in y, # in z)
	//int m_nGPTot;         //total# of gp

	int m_nTetNodNum;
	int m_nTetNum;
	int* m_npTetSeq;
	cor* m_spTetCor; // info about tetradron, which is used to judge gps' position

	gp* m_ptrGP;          // pointers to current gp 
	bgMesh* m_ptrMes;         // pointers to current mesh
	double* m_dpTrans;        // ultimate transfer function
	double* m_dpTrans_arranged;
	int m_nNodNum;          //# nodals
};

#endif 
