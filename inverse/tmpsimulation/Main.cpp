#include "inclu.h"
#include "matEg.h"
#include "Matrix.h"
#include "external.h"
#include "TMPSimulation.h"

matEg engine(MAX_BUF);
int ID;
int YDIM;
double DMAX;
int NUM_INF;
double BG_INTERVAL;
int IS_ISO;
int IS_HOMO;
char FHN_TYPE;
char PAR_TYPE;
int SCHEME;
char PATH_ROOT[MAX_FILENAME];



int main()
{
	readPara();
	
	
	int c;
	cout<<"Select function: 1, construct M,K; 2, TMP simulation"<<endl;
	cin>>c;
	
	if (c==1)
	{
		cout<<"Start constructing M,K"<<endl;
		Heart heart(PATH_ROOT);
		
		cout<<"Start assembling MK"<<endl;
		
		heart.assTrans();
		
		cout<<"M, K construction completed ^_____^"<<endl;
		heart.writeTrans();
		cout<<"M, K written ^_____^"<<endl;
	}
	else if (c == 2)
	{
		cout<<"Open matlab engine!"<<endl;
		engine.Open(NULL);
		cout<<"Constructing Simulator"<<endl;
		Simulator simul(PATH_ROOT);
		cout<<"Start simulationn"<<endl; 
		//clock_t s, f;
		time_t s = time(NULL);
		simul.staPro();
    	time_t f = time(NULL);
    	double d = difftime(f, s);
    	cout<<"time for A-P simulation "<<d<<endl;
    	
		cout<<"TMP simulation completed ^_______^"<<endl;
		simul.write(PATH_ROOT);
		cout<<"Results written ^___^"<<endl;
		engine.Close();
    	cout<<"Close matlab engine!"<<endl;		

	}
	else
		cout<<"Wrong command"<<endl;
	
	return EXIT_SUCCESS;
	
/*	char* file = new char[MAX_FILENAME];
	char* testPath = "/Users/maomaowlw/Research/Mywork/Workspace/TMPSimulation/Test/";
	
	//check mfree nodes - infR
	strcpy(file,testPath);
	strcat(file,"InfR.bin");
	heart.checkInfR(file);
	
	heart.checkInfMfree(testPath);
	
	//check gpts - cor    
	heart.checkGP(testPath);	
	
	delete []file;
*/	
}
