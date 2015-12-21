/*
 * external.cpp
 *
 *  Created on: Jul 3, 2012
 *      Author: jxx2144
 */

#include "inclu.h"
#include "external.h"
#include "Matrix.h"
#include "Heart.h"

//**********************************
//inVol
//*********************************
int inVol(tPointd a, tPointd b, tPointd c, tPointd d)
{
   double vol;
   double bxdx, bydy, bzdz, cxdx, cydy, czdz;

   bxdx=b[0]-d[0];
   bydy=b[1]-d[1];
   bzdz=b[2]-d[2];
   cxdx=c[0]-d[0];
   cydy=c[1]-d[1];
   czdz=c[2]-d[2];
   vol =   (a[2]-d[2]) * (bxdx*cydy - bydy*cxdx)
         + (a[1]-d[1]) * (bzdz*cxdx - bxdx*czdz)
         + (a[0]-d[0]) * (bydy*czdz - bzdz*cydy);

   //cout << "vol: " << vol << endl;

   /* The volume should be an integer. */
   if      ( vol > 0 )   return  1;
   else if ( vol < 0 )  return -1;
   else                    return  0;
}

//**************************
//in/outside a bunch of tetras
//***************************
bool outTetra(tPointd pt, double* cor, int* tet, int num)
{
	//loop through all tetras, to find whether at least one include this gp
	int pos = 0;
	for (int i=0; i<num; i++)
	{
		int n0 = tet[pos++];
		int n1 = tet[pos++];
		int n2 = tet[pos++];
		int n3 = tet[pos++];

		int pos1 = 3*(n0-1);
		tPointd c0,c1,c2,c3;
		c0[0]= cor[pos1];
		c0[1]= cor[pos1+1];
		c0[2]= cor[pos1+2];
		pos1 = 3*(n1-1);
		c1[0] = cor[pos1];
		c1[1] = cor[pos1+1];
		c1[2] = cor[pos1+2];
		pos1 = 3*(n2-1);
		c2[0] = cor[pos1];
		c2[1] = cor[pos1+1];
		c2[2] = cor[pos1+2];
		pos1 = 3*(n3-1);
		c3[0] = cor[pos1];
		c3[1] = cor[pos1+1];
		c3[2] = cor[pos1+2];

	    int V0 = inVol(c0, c1, c2, c3);
	    int V1 = inVol(pt, c1, c2, c3);
	    int V2 = inVol(c0, pt, c2, c3);
	    int V3 = inVol(c0, c1, pt, c3);
	    int V4 = inVol(c0, c1, c2, pt);

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

//*****************************************
//generate mfree particles, using class of "2DHeart"
//****************************************
void genMFreeNodes(char* file_in,char* file_out)
{
//	visualizer visu(file_in);
//	visu.saveImage(file_image,1);

	cout<<"Read in tetra file"<<endl;
	Heart2D heart(file_in);

	cout<<"Generate meshfree nodes"<<endl;
	heart.genMFreeNodes(file_out);

	cout<<"write results"<<endl;
	//double* mfree = heart.getResult();
	//int num = heart.getNumNodes();
//	visu.display2D(mfree,num);
//	visu.saveImage(file_image,2);
	wait();
}

void wait()
{
	char c;
	cout<<endl<<"Press any key to close the application"<<endl;
	cin>>c;

	return;
}
