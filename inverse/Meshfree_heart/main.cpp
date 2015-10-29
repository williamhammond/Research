/*
 * main.cpp
 *
 *  Created on: Jul 3, 2012
 *      Author: jxx2144
 */




#include "inclu.h"
#include "external.h"
#include "matEg.h"

matEg* engine = new matEg(MAX_BUFFER);
int main()
{
	cout<<"Meshfree particles generation"<<endl;

	char* path_in = "/Research/Meshfree/JHUPigData/Ref2/";
			//"/Users/maomaowlw/Research/Mywork/Data/Geometry/ECGSIM/NormalMaleYoung2/Ref/";
		//			"Auckland/Auckland_f02_2_5/Surfaces/";
    char* path_out = "/Research/Meshfree/JHUPigData/Final/";
    		//"/Users/maomaowlw/Research/Mywork/Data/Geometry/ECGSIM/NormalMaleYoung2/Final/";
		//			"PhysioNet/case4/Processed/Results/Ref/";
	char* filename_in = "heart";
	char* filename_out = "heart";
	//char* image = "Figures/heart3D_initial";
	char* file_in = new char[MAX_FILENAME];
	char* file_out = new char[MAX_FILENAME];
	//char* file_image = new char[MAX_FILENAME];


	strcpy(file_in,path_in);
	strcat(file_in,filename_in);
	strcpy(file_out,path_out);
	strcat(file_out,filename_out);
	//strcpy(file_image,path_out);
	//strcat(file_image,image);

    genMFreeNodes(file_in,file_out);
	delete []file_in;
	delete []file_out;

	//cout<<"Finish the Meshfree particles generation"<<endl;
	//delete []file_image;

}
