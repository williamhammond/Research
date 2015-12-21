/*
 * inclu.h
 *
 *  Created on: Jul 3, 2012
 *      Author: jxx2144
 */

#ifndef INCLU_H_
#define INCLU_H_

//************************************************************************
//  Includes
//************************************************************************
//#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//include <memory.h>
//#include <malloc.h>
#include <time.h>
#include<new>
#include <iostream>
#include <fstream>
//#include <tiffio.h>
using std::cerr;
using std::cout;
using std::cin;
using std::endl;


//************************************************************************
//  structures definition (to adopt Ken's function
//************************************************************************
//********************************
#define density 11.3
//#define dm 1.5 // influence domain scale
#define dm 3 // influence domain scale

#define PI 3.141592654
#define density 11.3
typedef int tPointi[3];
typedef double tPointd[3];

const int MAX_SLICE = 30;    //any value larger than the number of slices
const int MAX_NUM = 10000;   //any value larger than the number of nodes
const int MAX_BUFFER = 5e5;
const int MAX_K = 100000;
const int MAX_PM =200;         //max # of pacemakers
const int MIN_PMR = 2;       //min # of pacemakers on RV
const int MAX_PMR = 4;
const int MAX_PS = 100;        //max # of points per slice
const int MAX_MFREE = 10000;
const int MAX_FILENAME = 500;
const int MAX_ID = 5;
const int MAX_SUR = 10000;
const int MAX_SUPPORT = 400;
const int MAX_NB = 10;    //max # of neighbors of each surface node ( neighbors share the same line of a tri)
const int MAX_NB_LV = 8;
const double MAX_DIF_NOR = cos(15*PI/180);
const double MAX_DIS_LP = 0.1;      //maximum distance ratio between mfree nodes and P1P2
const int STEP_LOR = 3;    //sometimes need repeat several times to refine LV/RV (trial -- error)
const double MAX_R = 1000;
const int NB_MESH = 12;
const int DIM = 3;
const int NB_INFLU = 20;//16;    //whereever change this value, need to regenerate .domain file, so need to re-run gp generation
const int DM = 3;//1.5;//
const double SCALE = 100;
const double damping1 = 0.1;
const double damping2 = 0.1;
//const int ID_APEX = 8;
//const double u_mean_thres = 1e2;
const double delta_t =1e5;
//const bool enforce_obs = true;
const int IS_ANISO = 0;
const int NUM_SEG = 17;
const int INT_SEG = 2;
const int EXCITE = 1;
const double SLICE_INT = 10;
const int FIB_READY = 0;
const int SEG_READY = 1;
const double THRESH_RATIO_ENDO2EPI = 1.0/2.0;//1.0/1.5;//
const double THRESH_RATIO_L2RENDO = 1.0/2.0;
//for surface smoothing
const double LAMDA1 = 0.33;// 0.33;
const double LAMDA2 = -0.35;//0;//
const int STEP_SMOOTH = 100;
//for rotate the registered_fiber back to original system
const double RX1 = -20;//0; //-135;
const double RZ = 45;//0; //-45;
const double RX2 = 135;//0; //20;
//***for fibMapping*******************
const double POLY_A = 0.8187;//0.8132;//
const double POLY_B = 0.2042;//0.2119;//

//**for vtk************************
static tPointd COLOR_SUR_TARGET = {0,1,0};
static tPointd COLOR_SUR_REF = {0,0,1};
static tPointd COLOR_MFREE_TARGET = {1,0,0};
static tPointd COLOR_MFREE_REF = {1,1,0};
const double SIZE_POINT_MFREE = 3;
const double SIZE_POINT_SUR = 5;
const int ICP_TYPE = 2;
const int ICP_STEP = 10;

//***for writing out testing data********
extern char* PATH_TEST;
extern char* PATH_TEST_REF;


//*******************************
//segment

enum Seg{BaseAnte,BaseAnteSept,BaseInfeSept,BaseInfe,BaseInfeLate,BaseAnteLate,
         MidAnte,MidAnteSept,MidInfeSept,MidInfe,MidInfeLate,MidAnteLate,
		 ApicAnte,ApicSept,ApicInfe,ApicLate,
		 Apex};

const Seg SEG_ECTOPIC = MidAnte;//Apex;//MidAnteLate;//MidAnte;
class Matrix;


//*******************************
// gauss_pts
static double ri[] = {
0.0,
0.577350269189626, -0.577350269189626,
0.774596669241483, -0.774596669241483, 0,
0.861136311594053, -0.861136311594053, 0.339981043584856, -0.339981043584856,
0.906179845938664, -0.906179845938664, 0.538469310105683, -0.538469310105683, 0.0,
0.932469514203152, -0.932469514203152, 0.661209386466265, -0.661209386466265, 0.238619186083197, -0.238619186083197,
};

// gauss weight
static double alphai[] = {
2.0,
1.0, 1.0,
0.555555555555556, 0.555555555555556, 0.888888888888889,
0.347854845137454, 0.347854845137454, 0.652145154862546, 0.652145154862546,
0.236926885056189, 0.236926885056189, 0.478628670499366, 0.478628670499366, 0.568888888888889,
0.171324492379170, 0.171324492379170, 0.360761573048139, 0.360761573048139, 0.467913934572691, 0.467913934572691,
};

//*****************************
//   The definition of the Tree root node
//
typedef struct tree{
  struct nd *rootptr;
  int dims;
} Tree;
//*****************************
typedef struct {
	double w, w_dx, w_dy, w_dz;
} Weight;
//******************************
typedef struct {
	tPointd coor; // coordinates (x,y,z)
	tPointd fib; // fiber orientation (x,y,z)
//	tPointd normal; // normal
	double domain; // influence domain of a node
	bool is_sur;
	bool is_epi;
	bool is_lendo;
	bool is_rendo;
//	bool is_sep;
//	bool is_base;
//	bool is_apex;
//	bool is_fix;
//	bool is_obs;
	int is_excite;
	int segment;
	int slice;
} Node;    //general, suitable for surface nodes

//*******************************
typedef struct {
	tPointd coor;
	tPointd fib;
	int segment;
	int is_excite;
	bool is_lv;
	bool is_rv;
	//transmural position of the node
	int cell_trans;     //0=epi, 1 = endo; 2 = M  (small to large APD)
	double ratio_trans;
}Node_mfree;  //particularly for mfree nodes, used in simulation

//*************************************
/*
typedef struct {
	double* fib_ang_h, *fib_ang_v; // fiber angle (polar, azimuthal)
} Fiber;
*/

//*************************************
typedef struct {
	int nb_tri, nb_tetra, nb_epi_tri, nb_endo_tri;
	int *tri; // nodes of triangles
	int *tetra; // nodes of tetrahedra
//	int *epi_tri;
//	int *endo_tri;
//	double *tri_norm; // for tri
//	double *area; // for tri
//	bool *is_epi; // for tri
//	bool *is_endo; // for tri
} Elements;

//*************************************
typedef struct {
	int nb_epi, nb_lendo,nb_rendo,nb_sur; // nb_ebound1, nb_ebound2, nb_ebound, nb_fbound,
//	int *ebound1; // indices to the essential boundary nodes
//	int *ebound2;
//	int *fbound;
	int *epi;
	int *lendo;
	int *rendo;
	int *sur; // nodes on the surface
	double *surf; //*ebound1f, *ebound2f, *fboundf, *surf; // storing the values to be enforced
} Bound;
//*************************************
/*
typedef struct {
	int nb_obs_u, nb_obs_v;
	int *obs_u, *obs_v;
	double *obs_uf, *obs_vf;
} Observe;
*/

//will not allocate/free memory, ptrs are pointed to the memory

//*************************************
typedef struct {
	double a0, a1, a2, a3, a4, a5, a6, a7;
} NewM;
//*************************************
typedef struct {
	int nb_gpts_dim;
	int nb_mesh_x, nb_mesh_y, nb_mesh_z;
	double len; // each mesh is a square
	double xl, xu, yl, yu, zl, zu; // lower and upper bounds in different directions
} Bg_mesh;
//*************************************
/*
typedef struct {
	double Exx, Eyy, Ezz, Exy, Exz, Eyz;
} Strain;
*/
//*************************************

typedef struct {
	int id_epi, id_lendo, id_rendo;
	double r_epi, r_lendo, r_rendo;
	double dr,tr;   //tr = P1-mfreePts, dr = |P1P2|
} closestSur;





#endif /* INCLU_H_ */
