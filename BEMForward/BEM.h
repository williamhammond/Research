//BEM.h ---- declaration of bem class

#ifndef _BEM__H_
#define _BEM__H_
//***********************************************************
//basic class: BEM
//***********************************************************
class BEM
{
	//functions
public:
	// ctors/dtor
	BEM(char** filNamBEM, char** filNamMFree, CellType type);
	BEM();
	virtual ~BEM();

	// interface with users, include the all process of BEMForward
	virtual void BEMCompute(char** filOutput) = 0;

	double* getBEMK();
	double* getBEMP();
	int getNumSur();
	int getNumVol();

	//declaration of friends
	friend class BEMCell;
	friend class Tri3Cell;
	friend class MFree;

protected:
	virtual void intPK_first(int, int) = 0;
	virtual int isOutward() = 0;
	void intPK(int, int);
	virtual void eleNor(int)=0;
	int isBel(int nNodeSeq, int nEleSeq);
	void init(char** filNamMFree);
	bool isForbidden(int nodSeq);      //nodSeq from 1...N
	void rearrangePK();
	bool isDiscarded(int seq);      //seq is the index of array, so start from 0

	//attributes

	// the following ones construct the BEM object, get value through initializaton
	// for BEM	
protected:
	int m_nEleNum;
	int m_nNodNum;
	int* m_npNodSeq;
	double* m_dpNodCor;
//	double* m_dpNodCon;
	double* m_dpEleDer;
	double* m_dpEleAre;
	double* m_dpRefCor; //cor of inner origin, to test whether the normal computed is outward or inward
	cor varCor;
	double* m_dpGps;          //record down all gps.
	int m_nGPNum;
	int* m_npForbidden;       //nodes not used in elements
	int m_nForNum;               //# nodes not used
	int* m_npDiscarded;
	int m_nDiscarded;     //position of elements to be discarded in P,K matrices
	bool m_bisFirst;

	// the following ones are what needed in computing
	//for BEMCell
	double* m_dpK;
	double* m_dpP;
	double* m_dpK_arranged;
	double* m_dpP_arranged;
	Tri3Cell* m_ptrBEMEle;
	MFree* m_ptrMFree;
};
//***********************************************************
//derivative class I: BEM with MFree
//***********************************************************
class BEM_source: public BEM
{
	//functions
public:
	// ctors/dtor
	BEM_source(char** filNamBEM, char** filNamMFree, CellType type); //read in files for boundary and source, then
	//call BEM(file) and MFree(file) respectively
	virtual ~BEM_source();

	// interface with users, include the all process of BEMForward
	virtual void BEMCompute(char** filOutput);
	double* getBEMKa();
	double* getMFree();
private:
	virtual void intPK_first(int, int);
	virtual int isOutward();
	virtual void eleNor(int);
	virtual void rearrangePK();
	//declaration of friends
	friend class BEMCell;
	friend class Tri3Cell;
	friend class MFree;

	//attributes
private:
//	int* m_npFlag;                    // used to dicate whether this triangle has been picked up
	double* m_dpKa;
	double* m_dpKa_arranged;
	double* m_dpTetCor;
	int* m_npTet;          //tetra geo, used to detect in/out-ward normal
	int m_nTetNodNum;
	int m_nTetNum;
};

//*************************************************************
//derivative class II: BEM without source terms
//*************************************************************
class BEM_nosource: public BEM
{
	//functions
public:
	// ctors/dtor
	BEM_nosource(char** filNamBEM, char** filNamMFree, CellType type); //read in BEMfile(including two surfaces), then call BEM(double*,int*) to initialize
	virtual ~BEM_nosource();

	// interface with users, include the all process of BEMForward
	virtual void BEMCompute(char** filOutput);
private:
	virtual void intPK_first(int, int);
	virtual int isOutward();
	virtual void eleNor(int);
//	mwArray& getBEMKInner();
//    mwArray& getBEMK();
//	mwArray& getBEMKInner();
	//declaration of friends
	friend class BEMCell;
	friend class Tri3Cell;
private:
	int m_nNodNumInner;
	int m_nNodNumOuter;
	int m_nEleNumInner;
	int m_nEleNumOuter;
	int m_nForNumInner;
	int m_nForNumOuter;
//	mwArray m_dpKInner;
//	mwArray m_dpKOuter;
//	mwArray m_dpPInner;
};

#endif
