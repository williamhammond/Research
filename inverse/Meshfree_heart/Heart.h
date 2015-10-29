/*
 * Heart.h
 *
 *  Created on: Jul 3, 2012
 *      Author: jxx2144
 */

#ifndef HEART_H_
#define HEART_H_

#include "Matrix.h"
class Point;
//class Matrix;
class visualizer;
class localizer;
class vecRot;



//*************************************************************
// 2DHeart -- heart surface geometry, used to generate meshfree nodes
//            with different density
//*************************************************************
class Heart2D{

	friend class Heart;
	//*****************************
public:
	Heart2D(char* file_in);
	~Heart2D();

	void genMFreeNodes(char* file_out);
	double* getResult();
	int getNumNodes();

private:
	void bgSampling();
	bool outTetra(tPointd);
	void setDensity (int int_x = 8, int int_y = 8, int int_z = 10);

	//members********************
private:
	double* mdp_node_sur;
	int mn_num_node_sur;
	int* mnp_tetra;
	int mn_num_tet;
	double* mdp_node_mfree;            //final mfree nodes after getting rid of out-boundary nodes.
	int mn_num_node_mfree;

	tPointd mtp_bg_l, mtp_bg_u, mtp_bg_d;    //lower, upper and distance
	tPointd mtp_den;             // interval length in each direction
	int mnp_len[3];
	int mn_scheme;                   // how to generate mfree nodes.

};


#endif /* HEART_H_ */
