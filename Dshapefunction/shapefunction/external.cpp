/*
 * external.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: jxx2144
 */

//external.cpp
//definition of external funtions
#include "inclu.h"
//#include"BEMCell.h"
#include "external.h"
//#include "BEM.h"
#include "matEg.h"
#include "Matrix.h"
#include "gp.h"

//*********************************************************
//external variables
//*********************************************************
//char* FILE_REF_COR = "/Users/maomaowlw/Research/Mywork/Data/Geometry/Auckland/Auckland_f01/Auckland_f01.cor";
		//"";PhysioNet/case2/Processed/Results/Ref/heart.cor
//char*  FILE_REF_FIB = "/Users/maomaowlw/Research/Mywork/Data/Geometry/Auckland/Auckland_f01/Auckland_f01.fib";
		//"";PhysioNet/case2/Processed/Results/Ref/heart.fi

//***********************************************************
//tools
//***********************************************************
void addPath(char** fullname,char* path, char** name, int num)
{
	for(int i=0; i<num; i++)
	{
		if(name[i])
		{
			fullname[i] = new char[MAX_FILENAME];
    		strcpy(fullname[i],path);
    		strcat(fullname[i],name[i]);
		}
		else
			fullname[i] = NULL;
	}
}

void wait()
{
	cout<<"Press any key to continue"<<endl;
	char c;
	cin>>c;
	exit(0);
}

void reSet(double* data, int len)
{
	for(int i=0; i<len; i++)
		data[i] = 0;
}
//-*******************************************************
//inVol
//called by "isout"
//*******************************************************
int inVol(cor a, cor b, cor c, cor d)
{
   double vol;
   double bxdx, bydy, bzdz, cxdx, cydy, czdz;

   bxdx=b.x-d.x;
   bydy=b.y-d.y;
   bzdz=b.z-d.z;
   cxdx=c.x-d.x;
   cydy=c.y-d.y;
   czdz=c.z-d.z;
   vol =   (a.z-d.z) * (bxdx*cydy - bydy*cxdx)
         + (a.y-d.y) * (bzdz*cxdx - bxdx*czdz)
         + (a.x-d.x) * (bydy*czdz - bzdz*cydy);

   //cout << "vol: " << vol << endl;

   /* The volume should be an integer. */
   if      ( vol > 0 )   return  1;
   else if ( vol < 0 )  return -1;
   else                    return  0;
}

//*******************************************************
//isout
//judge whether gp is out of volume
//*******************************************************
bool isoutTetra(cor gpCor,double* node, int* tet, int num_tet)
{
	//loop through all tetras, to find whether at least one include this gp
	int pos = 0;
	int pos1 = 0;
	for (int i=0; i<num_tet; i++)
	{
		int n0 = tet[pos++];
		int n1 = tet[pos++];
		int n2 = tet[pos++];
		int n3 = tet[pos++];

		cor c0;
		pos1 = 3*(n0-1);
		c0.x= node[pos1++];
		c0.y = node[pos1++];
		c0.z = node[pos1];
		cor c1;
		pos1 = 3*(n1-1);
		c1.x = node[pos1++];
		c1.y = node[pos1++];
		c1.z = node[pos1];
		cor c2;
		pos1 = 3*(n2-1);
		c2.x = node[pos1++];
		c2.y = node[pos1++];
		c2.z = node[pos1];
		cor c3;
		pos1 = 3*(n3-1);
		c3.x = node[pos1++];
		c3.y = node[pos1++];
		c3.z = node[pos1];


	   int V0 = inVol(c0, c1, c2, c3);
	   int V1 = inVol(gpCor, c1, c2, c3);
	   int V2 = inVol(c0, gpCor, c2, c3);
	   int V3 = inVol(c0, c1, gpCor, c3);
	   int V4 = inVol(c0, c1, c2, gpCor);

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

//*********************************************************
//reload function to compute cos(theta) between two vetors)
//*********************************************************
//ordinary
double dotProduct(double x,double y, double z, double x1, double y1, double z1)
{
	return (x*x1 + y*y1 + z*z1)/(sqrt(x*x + y*y + z*z)*sqrt(x1*x1+y1*y1+z1*z1));
}

//with x,y or z axis, simplified
double dotProduct(char axis, double x, double y, double z)
{
	switch (axis)
	{
	case 'x':
		return x/sqrt(x*x + y*y + z*z);
	case 'y':
		return y/sqrt(x*x + y*y + z*z);
	case 'z':
		return z/sqrt(x*x + y*y + z*z);
	default:
		cout<<"error: nonexistent cartesian basis!"<<endl;
		exit(0);
	}
}

//**********************************************************
//Poly3D family for p, pdx...
//**********************************************************
void Poly3D(Matrix& p, double x, double y, double z)
{
	p.setData(1.0,0,0);
	p.setData(x,1,0);
	p.setData(y,2,0);
	p.setData(z,3,0);
}

void Poly3D_dx(Matrix& p, double x, double y, double z)
{
	p.setData(0.0,0,0);
	p.setData(1.0,1,0);
	p.setData(0.0,2,0);
	p.setData(0.0,3,0);
}

void Poly3D_dy(Matrix& p, double x, double y, double z)
{
	p.setData(0.0,0,0);
	p.setData(0.0,1,0);
	p.setData(1.0,2,0);
	p.setData(0.0,3,0);
}

void Poly3D_dz(Matrix& p, double x, double y, double z)
{
	p.setData(0.0,0,0);
	p.setData(0.0,1,0);
	p.setData(0.0,2,0);
	p.setData(1.0,3,0);
}



//***************************************
//check whether a matrix has nan
//****************************************
int checkNan(Matrix& p)
{
	double* ptr = p.getDataPtr();
	int num = p.getCol()*p.getRow();
	for(int i=0; i<num; i++)
	{
		if (isnan(ptr[i]))
			return 1;
	}
	return 0;
}


void computeDeltaFECG()
{
	    int NUM_NODE = 1401;
		Matrix finalshape;
		finalshape.setData(0.0,3*NUM_NODE,NUM_NODE);
		MFree* mfree;
		char* file = new char[MAX_FILENAME];
		char* fileout = new char[MAX_FILENAME];
		char* copath_out = new char[MAX_FILENAME];
		char id[5];
		sprintf(id,"%d",NUM_NODE);
		//strcpy(file, "/////");

	   char* path ="/Research/Meshfree/AW/Final/1401/";
			//"/Research/delta_newH/";
			//"/users/axr8834/desktop/cube729/";
	   strcpy(file,path);
	   strcat(file,"heart_mfree_1401.cor");
	   copath_out = "/Research/Meshfree/AW/Simulation/1401/Input/";

		mfree = new MFree(file);
	//	mfree->m_nGPNum = m_nMFNum;
		mfree->assTrans(file);

		for (int j = 0;j<NUM_NODE;j++){

			finalshape.setSub(mfree->m_mpGPshaD[j],3*j,3*j+2,0,mfree->m_npGPNBNum[j]-1);
		    }
		cout<<"time for shape functions"<<endl;


		char** filNamDelta = new char*[3];
		char** namdelta = new char*[3];
		namdelta[0] = "delta.bin"; //"heart_mfnode_azar_4fibermapping.cor";//"cube1.cor";//"sphCor_mfree_1388.bin";//"auckland_1720.cor";//
		namdelta[1] = "index.bin"; //"heart_mfnode_azar_4fibermapping.fib";//"cube1.fib";//"sphFib_mfree_1388.bin";// "auckland_1720.fib"; //
		namdelta[2] = "surnum.bin";//"heart_mfnode_azar_4fibermapping_tet.cor";//"cube1.tet.cor";//"sphCor_sur_inner_1136.bin";// "auckland_sur1254.cor";//"human_sur1911.cor";

		addPath(filNamDelta,copath_out,namdelta,3);

		writeFile(filNamDelta[0],finalshape.getDataPtr(),sizeof(double),3*NUM_NODE*NUM_NODE,"wb");
		writeFile(filNamDelta[1],mfree->m_npGPNB,sizeof(int),NUM_NODE*MAX_INF,"wb");
		writeFile(filNamDelta[2],mfree->m_npGPNBNum,sizeof(int),NUM_NODE,"wb");
		delete mfree;
		//delete []file;
}



