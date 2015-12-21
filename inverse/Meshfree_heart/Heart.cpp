/*
 * Heart.cpp
 *
 *  Created on: Jul 3, 2012
 *      Author: jxx2144
 */


#include "inclu.h"
#include "external.h"
#include "Heart.h"
#include "matEg.h"
#include "Matrix.h"
//extern matEg* engine;

//******************************************************************
//Class: 2DHeart
// Used to generate mfree nodes given heart surface geometry
//******************************************************************
//***********************
// ctor(s)/dtor
//**********************
Heart2D::Heart2D(char* file_in)
{
	char* file = new char[MAX_FILENAME];
	strcpy(file,file_in);
	strcat(file,"_tet.cor");
	readFile(file,&mdp_node_sur,sizeof(double),3,mn_num_node_sur);
	strcpy(file,file_in);
	strcat(file,".tet");
	readFile(file,&mnp_tetra,sizeof(int),4,mn_num_tet);

	mdp_node_mfree = NULL;
	mn_num_node_mfree = 0;
}

Heart2D::~Heart2D()
{
	if(mdp_node_sur)
		delete []mdp_node_sur;
	if(mdp_node_mfree)
		delete []mdp_node_mfree;
	if(mnp_tetra)
		delete []mnp_tetra;
}

//*********************
//set how many intervals should be in each direction
//********************
void Heart2D::setDensity(int int_x, int int_y, int int_z)
{
	//first compute the dx,dy,dz

	mtp_bg_l[0] = mdp_node_sur[0];
	mtp_bg_l[1] = mdp_node_sur[1];
	mtp_bg_l[2] = mdp_node_sur[2];
	mtp_bg_u[0] = mdp_node_sur[0];
	mtp_bg_u[1] = mdp_node_sur[1];
	mtp_bg_u[2] = mdp_node_sur[2];
	//loop through every tetra nodes to determine the lower bound and upper bound of bg mesh
	for(int i=0; i< mn_num_node_sur; i++)
	{
		if(mtp_bg_l[0] >mdp_node_sur[3*i])
			mtp_bg_l[0] = mdp_node_sur[3*i];
		if(mtp_bg_l[1] > mdp_node_sur[3*i+1])
			mtp_bg_l[1]  = mdp_node_sur[3*i+1];
		if(mtp_bg_l[2] > mdp_node_sur[3*i+2])
		    mtp_bg_l[2]  = mdp_node_sur[3*i+2];
		if(mtp_bg_u[0] < mdp_node_sur[3*i])
			mtp_bg_u[0]= mdp_node_sur[3*i];
		if(mtp_bg_u[1] < mdp_node_sur[3*i+1])
			mtp_bg_u[1] = mdp_node_sur[3*i+1];
		if(mtp_bg_u[2] <mdp_node_sur[3*i+2])
			mtp_bg_u[2] = mdp_node_sur[3*i+2];
	}

	//get the directional #  -- round off into higher integral
	mtp_bg_d[0] = mtp_bg_u[0] - mtp_bg_l[0];
	mtp_bg_d[1] = mtp_bg_u[1] - mtp_bg_l[1];
	mtp_bg_d[2] = mtp_bg_u[2] - mtp_bg_l[2];

    //first display the dimension of bgmesh
    cout<<"*******************************************************"<<endl;
	cout<<"bg mesh dimension (x, y, z): "<<mtp_bg_d[0]<<"  "<<mtp_bg_d[1]<<"  "<<mtp_bg_d[2]<<endl;
	cout<<"******************************************************"<<endl;
	cout<<"Need adjust density (number of elements in each direction) ? (y/n)"<<endl;
	char flag;
	cin>>flag;
	if (flag == 'y' || flag == 'Y')
	{
		cout<<"x:  ";
		cin>>int_x;
	    cout<<endl;
	    cout<<"y:  ";
	    cin>>int_y;
	    cout<<endl;
	    cout<<"z:  ";
	    cin>>int_z;
	}
	mnp_len[0] = int_x;
	mnp_len[1] = int_y;
    mnp_len[2] = int_z;

	//compute density
	mtp_den[0] = ceil(mtp_bg_d[0]/ int_x);
	mtp_den[1] = ceil(mtp_bg_d[1]/ int_y);
	mtp_den[2] = ceil(mtp_bg_d[2]/ int_z);
	// adjust the boundary cor to make integral# mesh
	double d = mtp_den[0]*int_x - mtp_bg_d[0];
	mtp_bg_l[0] -= d/2.0;
	mtp_bg_u[0] += d/2.0;
	d = mtp_den[1]*int_y - mtp_bg_d[1];
	mtp_bg_l[1] -= d/2.0;
	mtp_bg_u[1] += d/2.0;
	d = mtp_den[2]*int_z - mtp_bg_d[2];
	mtp_bg_l[2] -= d/2.0;
	mtp_bg_u[2] += d/2.0;

	mdp_node_mfree = new double[int_x*int_y*int_z*DIM];

}

//**********************************
//test wheter a sample point is in the tetra
//***********************************
bool Heart2D::outTetra(tPointd cor)
{
	//loop through all tetras, to find whether at least one include this gp
	int pos = 0;
	for (int i=0; i<mn_num_tet; i++)
	{
		int n0 = mnp_tetra[pos++];
		int n1 = mnp_tetra[pos++];
		int n2 = mnp_tetra[pos++];
		int n3 = mnp_tetra[pos++];

		int pos1 = 3*(n0-1);
		tPointd c0,c1,c2,c3;
		c0[0]= mdp_node_sur[pos1];
		c0[1]= mdp_node_sur[pos1+1];
		c0[2]= mdp_node_sur[pos1+2];
		pos1 = 3*(n1-1);
		c1[0] = mdp_node_sur[pos1];
		c1[1] = mdp_node_sur[pos1+1];
		c1[2] = mdp_node_sur[pos1+2];
		pos1 = 3*(n2-1);
		c2[0] = mdp_node_sur[pos1];
		c2[1] = mdp_node_sur[pos1+1];
		c2[2] = mdp_node_sur[pos1+2];
		pos1 = 3*(n3-1);
		c3[0] = mdp_node_sur[pos1];
		c3[1] = mdp_node_sur[pos1+1];
		c3[2] = mdp_node_sur[pos1+2];

	    int V0 = inVol(c0, c1, c2, c3);
	    int V1 = inVol(cor, c1, c2, c3);
	    int V2 = inVol(c0, cor, c2, c3);
	    int V3 = inVol(c0, c1, cor, c3);
	    int V4 = inVol(c0, c1, c2, cor);

	if( V0 == 0 )
	{
	//	cerr << "The tetrahedron is degenerate !!" << endl;
	//	exit(0);
	//	return true;
	}
	else
	{
		if( V0 > 0 )
		{
		 if( V1>=0 && V2>=0 && V3>=0 && V4>=0)
				return false; // true only if all are same sign as V0
		}	//else
			//	return true;
	   else
	   {
			if( V1<=0 && V2<=0 && V3<=0 && V4<=0)
				return false; // true only if all are same sign as V0
	   }	//	else
		//		return true;
	}
	}
	return true;
}

//*******************************
//generate meshfree particles with different scheme
//*******************************
void Heart2D::bgSampling()
{
	//choose which kind of scheme to be used
	//1, vertex of each cube
	//2, vertex plus central node
	cout<<"Pls select the scheme of sampling:"<<endl;
	cout<<"1, vertex only;          2, with central nodes"<<endl;
	cin>>mn_scheme;

	tPointd temp_cor;
	int pos = 0;
	for(int i=0; i<=mnp_len[2];i++)
	{
		for(int j=0; j<=mnp_len[1]; j++)
		{
			for(int k=0; k<=mnp_len[0]; k++)
			{
				temp_cor[0] = mtp_bg_l[0] + k*mtp_den[0];
				temp_cor[1] = mtp_bg_l[1] + j*mtp_den[1];
				temp_cor[2] = mtp_bg_l[2] + i*mtp_den[2];
				if(!outTetra(temp_cor))
				{
					pos = 3*mn_num_node_mfree;
					mdp_node_mfree[pos++] = temp_cor[0];
					mdp_node_mfree[pos++] = temp_cor[1];
					mdp_node_mfree[pos++] = temp_cor[2];
                    mn_num_node_mfree++;
				}
			}
		}
	}
	//additional nodes in scheme 2
	if (mn_scheme == 2)
	{
		for(int i=0; i<mnp_len[2]; i++)
		{
			for(int j=0; j<mnp_len[1];j++)
			{
				for(int k=0; k<mnp_len[0]; k++)
				{
					temp_cor[0] = mtp_bg_l[0] + (2*k+1)*mtp_den[0]/2.0;
					temp_cor[1] = mtp_bg_l[1] + (2*j+1)*mtp_den[1]/2.0;
					temp_cor[2] = mtp_bg_l[2] + (2*i+1)*mtp_den[2]/2.0;
					if(!outTetra(temp_cor))
					{
						pos = 3*mn_num_node_mfree;
						mdp_node_mfree[pos++] = temp_cor[0];
						mdp_node_mfree[pos++] = temp_cor[1];
						mdp_node_mfree[pos++] = temp_cor[2];
						mn_num_node_mfree++;
					}
				}
			}
		}
	}
	cout<<"Total # mfree nodes:  "<<mn_num_node_mfree<<endl;
}

//*********************************************
// public function to be called by other application
//******************************************
//************************
//generate a set of new mfree nodes
//************************
void Heart2D::genMFreeNodes(char* file_out)
{
	setDensity();
	bgSampling();

	char id[8];
	sprintf(id,"%d",mn_num_node_mfree);
	char* file = new char[MAX_FILENAME];
    strcpy(file,file_out);
	strcat(file,"_mfree_");
	strcat(file,id);
	strcat(file,".cor");
	writeFile(file,mdp_node_mfree,sizeof(double),3*mn_num_node_mfree,"wb");
	delete []file;
}



//**************************************
//obtain resulted mfree nodes
//**************************************
double* Heart2D::getResult()
{
	return mdp_node_mfree;
}

int Heart2D::getNumNodes()
{
	return mn_num_node_mfree;
}

