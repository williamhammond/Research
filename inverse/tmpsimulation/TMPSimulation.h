#ifndef TMPSIMULATION_H_
#define TMPSIMULATION_H_

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
	void setData(double, int, int);
public:
	Matrix d;
	Matrix dx;
	Matrix dy;
	Matrix dz;
};

//*******************************************************
//class: bgMesh
//the whole bgMesh info, recording the boundary of the mesh, also 
//used to generate gp in each gp mesh; judge whether outof original volume domain
//*******************************************************
class bgMesh {
public:
	bgMesh(double* cor, int num);
	~bgMesh();

public:
	void generateGPAll(double* cor, int* tet, int num, int num_tet);
private:
	void setMesSeq(int, int, int);
	void genGloGP(double* gpcor, double* gpwei, int l, int m, int n);
	void comSta(double* cor, int num); // 0, totally out; 1, partially out; else, in
	bool isOut(double* gpcor, double* cor, int* tet, int num_tet);

private:
	//member
	//global mesh info: boundary, number of mesh each direction idx of which mesh it is
	double m_dgxl;
	double m_dgxu;
	double m_dgyl;
	double m_dgyu;
	double m_dgzl;
	double m_dgzu;
	//the increase of the boundary on each direction, so that to fit in integer number of meshes.
	double m_dTolX;
	double m_dTolY;
	double m_dTolZ;
	//number of meshes in each direction
	int m_nNumX;
	int m_nNumY;
	int m_nNumZ;
	int m_ngi;
	int m_ngj;
	int m_ngk;

	//local mesh info: boundary of each mesh
	double m_dxl;
	double m_dxu;
	double m_dyl;
	double m_dyu;
	double m_dzl;
	double m_dzu;

	double m_dMesInt;
	int m_nSta;
	// shared info, only initialized once
	double* m_dpNatGPCor;
	double* m_dpNatGPWei;
	int m_nGPNum;         // # gp for each mesh
public:
	int m_nGPNumAll;
	double* m_dpGPCor;
	double* m_dpGPWei;
};

//*******************************************************
//class:nodal
//Abstract generic representation of a node, with coordinate, fiber
//*******************************************************
class node {
public:
	node();
	node(double*, double*, double D);
	~node();
	//re-intialization
protected:
	void reInit(double*, double*, double D);

	//member
protected:
	double* m_dpCor;
	double* m_dpFibOrWei;        //fib for mfree nodes, weight for gs poitns
	double m_dD;       // Difussion tensor value;
};

//*******************************************************
//class:meshfree point
//derivative of class node
//*******************************************************
class nodeMF: public node {
	friend class nodeGS;
	friend class Heart;
public:
	nodeMF();
	nodeMF(double*, double*, double, double*, int);
	~nodeMF();
	void reInit(double*, double*, double, double*, int);

private:
	void calculRInf(double* cor, int num);
private:
	int* m_npIInf;        //idx of nodes inside the influnce domain of this node
	double m_dRInf;         //radius of the 
	Matrix m_mp;
};

//*******************************************************
//class:gs point
//derivative of class node
//*******************************************************
class nodeGS: public node {
	friend class nodeMF;
	friend class Heart;
public:
	nodeGS();
	nodeGS(double*, double*, nodeMF*, int);
	~nodeGS();
	void reInit(double*, double*, nodeMF*, int);
	double getWei();
private:
	void initMatrixInf();
	void calInf(nodeMF* node, int num);
	void calSha(nodeMF* node);
	void calDB(nodeMF* node);
	void calRot();
	void calRotbyRotation();
	void calFib(nodeMF* node);
public:
	int m_nSur;     //# suppoting nodes
	int* m_npISur;   //id of supporting nodes
	Matrix m_mSha;
	Matrix m_mDB;
private:
	double* m_dpFib;
	double* m_dpWSur; //weight associated with supporting nodes, including w, dwdx, dwdy, dwdz
	MatrixGroup m_mP;
	MatrixGroup m_mA;
	MatrixGroup m_mB;
	MatrixGroup m_mInvA;
	Matrix m_mD;
	Matrix m_mQ;  //rotation matrix
};

//*****************************************************
//class: Heart
//including all necessary info of the heart: tet;meshfree nodes; gs nodes;
//*****************************************************
class Heart {
public:
	Heart(char* path);
	~Heart();
private:
	void readGeo();
	void readMF();
	void genBG();
	void genGS();
public:
	void assTrans();
	void writeTrans();

	void checkInfR(char* testPath);
	void checkGP(char* testPath);
	void checkInfMfree(char* file);

private:
	char* path;
	int m_nTetNum;      //# tetra
	int m_nTetVerNum;   //# vertice of tetra
	double* m_dpTetCor;
	int* m_npTet;
	int m_nMFNum;

	nodeMF* mfree;
	nodeGS* gspt;
	bgMesh* bg;

	double* m_dpM;
	double* m_dpK;
};

//********************************************************
//class: Simulator
//**********************************************************
class Simulator {
public:
	Simulator(char* path);
	~Simulator();
	void staPro();
	void write(char* path_out);
private:
	void genInitLV();
	void genInitRV();
	void genPacing();
	void processA(char lr);
	void removeSti(char);
	void addSti(char);

	void getCurrentSti(int i);
	void updateA(int i);
	void samPro_RK(int i, double dt);
	void fg(Matrix& u, Matrix& v);
	void setExcitation(char* path);
	void setPacing(char* in);

	double numDiff(double* head_u, double* head_t, char type);
	void iniAT(int& s_pos, int& f_pos, int& rs_pos, int& rf_pos, int& max_pos,
			double& max, double* u);
	void finPeak(int& peak_pos, int& re_pos, double* u);
	void u2t();

private:
	double* t;
	double tmax;
	double dt_simul;
	double dt_save;

	Matrix U;
	Matrix V;
//	Matrix Y;
	Matrix T;           // Trans Matrtix in TMP model
	Matrix H;           // Transfer Matrix in TMP-to-BSP model;
	Matrix Sti_all;
	Matrix Sti;
	int* sti_l;
	int* sti_r;
	int* sti_p;   //sti indexing whether l/r/p has stimulus at each time instant
	Matrix Par;
	Matrix state;
	Matrix mea;
	int step_simul;
	int step_save;
	double t_delay;
	int step_delay;
	int* exc;
	int num_exc;
	int num_exc_lv;
	int num_exc_rv;
	double duration_exc;

	double* focus;
	int num_focus;
	double duration_focus;
	double freq_focus;
	double start_focus;
	int startstep_focus;

	Matrix u0_RK;
	Matrix v0_RK;
	Matrix ui_RK;
	Matrix vi_RK;
	Matrix u0_fg;
	Matrix v0_fg;
	Matrix vec_one;

	double* parameter;

	double* AT;
	int* AT_idx;
	double* RT;
	int* RT_idx;
	double* APD;
};

#endif /*TMPSIMULATION_H_*/
