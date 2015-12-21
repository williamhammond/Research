// Filename: BEMCell.h
// declaration of classes used in this project

#ifndef __BEMCELL_H__
#define __BEMCELL_H__

//**********************************************************************************************************************
// declaration/definition of constant, enum,struct, and class
//**********************************************************************************************************************

//************************************************************
// represent the cell type used
enum CellType {Lin2,Tri3, Rec4, Tri6};
//************************************************************
//struct: represent different structures and parameters for different cell types
struct CellPara;
//************************************************************
// declaration of class type
class BEM;

//************************************************************************************************************************
//declaration of classes
//************************************************************************************************************************
//************************************************************
//Class name: BEMCell
//Class type: abstract basic class
//Function: provide an mutual interface for different kinds of cell types
//************************************************************

class BEMCell
{
	friend class BEM;
	friend class BEM_source;
	friend class BEM_nosource;
public:
	//ctor(s)/dtor
	BEMCell ( int nNodNum, int nDim);

	virtual ~BEMCell ();
	
	// member functions
    // set the sequence of processed node
    void setNodSeq ( int, const BEM& );
	// set the sequence of processd element
	void setEleSeq ( int, const BEM& );
	// compute single K,P, virtual function
	virtual void sinKP (BEM& bem, int isIn, int isFirst) = 0;
    // used to test whether the single integral on the BE 
	// is close to analytical solution
	virtual double test () = 0;

protected:
	//avoid = between different derivatives
	void operator= (const BEMCell&); 
	// numeric integral, virtual function
	virtual void numInt (BEM& bem, int isFirst) = 0;
	// numeric integral for sigular integral
	virtual void sinInt () = 0;
	// generate shape function
	virtual void gauQua_normal ( BEM& bem, int isFirst, double* data_K, double* data_P, double h,double s, double wr, double ws = 0) = 0;
	virtual void gauQua_sin (double* data_K, double* data_P, double h,double s, double wr, double ws = 0) = 0;
	// compute |G| --- include normal derivative, |G| = 2*area;
	virtual void norG () = 0;
	// generate guassian points and its weights
	virtual void genGau ( int, int) = 0; 

//	int isOutward(BEM& bem, char type);
   
	//attributes
protected:
	int m_nNodSeq;
	int m_nEleSeq;
	int m_nNodNum;
	int m_nDim;
	// value of |G|
	double m_nG;
	// norm derivatives of element
	double *m_dpNor;
	// nodal coordinates, 2D array
    double *m_dpNodCor;
	// nodel anisotropic conductance, for later extension
	double *m_dpNodCon;
	// point coordinate, 1D array
	double *m_dpPoiCor;
	// Single matrix for K&P ( actually vector)
	double *m_dpsinK;
	double *m_dpsinP;
	double* m_dpSha;       //shape function for each gp.
};


//************************************************************
//class name: Tri3Cell
//class type: derivative class of BEMCell
//function: triangular element with 3 nodes, linear element with natural coordinate L1,L2,L3=1-L1-L2
//************************************************************
class Tri3Cell:public BEMCell
{
	friend class BEM;
	friend class BEM_source;
	friend class BEM_nosource;
public:
	// ctors/dtor
	Tri3Cell(int nNodNum = 3, int nDim = 3, int nGauNum = GAUPOI,int nsinGauNum = SINGAUPOI);
	virtual ~Tri3Cell ();
  
	// functions
public:
	//virtual void getPara ( const FEM& );
    virtual void sinKP (BEM& bem, int isIn, int isFirst);
	
	virtual double test ();

private:
	// numeric integral, virtual function
	virtual void numInt (BEM& bem, int isFirst );
	void numInt();
	// sigular integral
	virtual void sinInt ();
	// compute |G| = 2A
	virtual void norG ();
	// compute area
	//generate Gaussian
	virtual void genGau (int,int);
	virtual void gauQua_normal ( BEM& bem, int isFirst, double* data_K, double* data_P, double h,double s, double wr, double ws = 0);
	virtual void gauQua_sin (double* data_K, double* data_P, double h,double s, double wr, double ws = 0);
	void gauQua_normal (double* data_K, double* data_P, double h,double s, double wr, double ws = 0);
    void gauQua_sur(double* data_S, double h,double s,double wr, double ws = 0);
	void surInt();
	void setEleSeq(int, const BEM&);
	//member
private:
	int m_nGauNum;
	int m_nSinGauNum;
	double *m_dpGauPoiCor;
	double *m_dpGauPoiWei;
	//for rectangle
	double *m_dpRecGauPoiCor;
	double *m_dpRecGauPoiWei;
	double *m_dpsinS;
};



#endif
