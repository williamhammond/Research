//external.cpp
//definition of external funtions
#include "inclu.h"
#include"BEMCell.h"
#include "external.h"
#include "MFree.h"
#include "BEM.h"
#include "matEg.h"
#include "Matrix.h"
//*********************************************************
//external variables
//*********************************************************
//char* FILE_REF_COR = "/Users/maomaowlw/Research/Mywork/Data/Geometry/Auckland/Auckland_f01/Auckland_f01.cor";
//"";PhysioNet/case2/Processed/Results/Ref/heart.cor
//char*  FILE_REF_FIB =  
//"/Users/maomaowlw/Research/Mywork/Data/Geometry/Auckland/Auckland_f01/Auckland_f01.fib";
//"";PhysioNet/case2/Processed/Results/Ref/heart.fi

//***********************************************************
//tools
//***********************************************************
void addPath(char** fullname, char* path, char** name, int num) {
	for (int i = 0; i < num; i++) {
		if (name[i]) {
			fullname[i] = new char[MAX_FILENAME];
			strcpy(fullname[i], path);
			strcat(fullname[i], name[i]);
		} else
			fullname[i] = NULL;
	}
}

void wriFile(char* file, cor* data, int num, char* mode) {
	FILE* fid = fopen(file, mode);
	double* d = new double[3];
	for (int i = 0; i < num; i++) {
		c2d(d, data[i]);
		fwrite(d, sizeof(double), 3, fid);
	}
	fclose(fid);
}

/*
 void wriFile(char* file, mwArray& mat,char* mode)
 {
 FILE* fid = fopen(file,mode);
 double* d = mxGetPr(mat.GetData());
 int row,col;
 row = size(&col,mat);
 fwrite(d,sizeof(double),row*col,fid);
 fclose(fid);
 }


 void test(char* filName,double* data,int len, char* mode)
 {
 char* testFile = new char[MAX_FILENAME];
 strcpy(testFile,"../BEMForward/test/");
 strcat(testFile,filName);
 wriFile(testFile,data,len,mode);
 }
 */

void wait() {
	cout << "Press any key to continue" << endl;
	char c;
	cin >> c;
	exit(0);
}

void reSet(double* data, int len) {
	for (int i = 0; i < len; i++)
		data[i] = 0;
}
//-*******************************************************
//inVol
//called by "isout" 
//*******************************************************
int inVol(cor a, cor b, cor c, cor d) {
	double vol;
	double bxdx, bydy, bzdz, cxdx, cydy, czdz;

	bxdx = b.x - d.x;
	bydy = b.y - d.y;
	bzdz = b.z - d.z;
	cxdx = c.x - d.x;
	cydy = c.y - d.y;
	czdz = c.z - d.z;
	vol = (a.z - d.z) * (bxdx * cydy - bydy * cxdx)
			+ (a.y - d.y) * (bzdz * cxdx - bxdx * czdz)
			+ (a.x - d.x) * (bydy * czdz - bzdz * cydy);

	//cout << "vol: " << vol << endl;

	/* The volume should be an integer. */
	if (vol > 0)
		return 1;
	else if (vol < 0)
		return -1;
	else
		return 0;
}

//*******************************************************
//isout
//judge whether gp is out of volume
//*******************************************************
bool isoutTetra(cor gpCor, double* node, int* tet, int num_tet) {
	//loop through all tetras, to find whether at least one include this gp
	int pos = 0;
	int pos1 = 0;
	for (int i = 0; i < num_tet; i++) {
		int n0 = tet[pos++];
		int n1 = tet[pos++];
		int n2 = tet[pos++];
		int n3 = tet[pos++];

		cor c0;
		pos1 = 3 * (n0 - 1);
		c0.x = node[pos1++];
		c0.y = node[pos1++];
		c0.z = node[pos1];
		cor c1;
		pos1 = 3 * (n1 - 1);
		c1.x = node[pos1++];
		c1.y = node[pos1++];
		c1.z = node[pos1];
		cor c2;
		pos1 = 3 * (n2 - 1);
		c2.x = node[pos1++];
		c2.y = node[pos1++];
		c2.z = node[pos1];
		cor c3;
		pos1 = 3 * (n3 - 1);
		c3.x = node[pos1++];
		c3.y = node[pos1++];
		c3.z = node[pos1];

		int V0 = inVol(c0, c1, c2, c3);
		int V1 = inVol(gpCor, c1, c2, c3);
		int V2 = inVol(c0, gpCor, c2, c3);
		int V3 = inVol(c0, c1, gpCor, c3);
		int V4 = inVol(c0, c1, c2, gpCor);

		if (V0 == 0) {
			//	cerr << "The tetrahedron is degenerate !!" << endl;
			//	exit(0);
			//	return true;
		} else {
			if (V0 > 0) {
				if (V1 >= 0 && V2 >= 0 && V3 >= 0 && V4 >= 0)
					return false; // true only if all are same sign as V0
			}	//else
				//	return true;
			else {
				if (V1 <= 0 && V2 <= 0 && V3 <= 0 && V4 <= 0)
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
double dotProduct(double x, double y, double z, double x1, double y1,
		double z1) {
	return (x * x1 + y * y1 + z * z1)
			/ (sqrt(x * x + y * y + z * z) * sqrt(x1 * x1 + y1 * y1 + z1 * z1));
}

//with x,y or z axis, simplified
double dotProduct(char axis, double x, double y, double z) {
	switch (axis) {
	case 'x':
		return x / sqrt(x * x + y * y + z * z);
	case 'y':
		return y / sqrt(x * x + y * y + z * z);
	case 'z':
		return z / sqrt(x * x + y * y + z * z);
	default:
		cout << "error: nonexistent cartesian basis!" << endl;
		exit(0);
	}
}

//**********************************************************
//Poly3D family for p, pdx...
//**********************************************************
void Poly3D(Matrix& p, double x, double y, double z) {
	p.setData(1.0, 0, 0);
	p.setData(x, 1, 0);
	p.setData(y, 2, 0);
	p.setData(z, 3, 0);
}

void Poly3D_dx(Matrix& p, double x, double y, double z) {
	p.setData(0.0, 0, 0);
	p.setData(1.0, 1, 0);
	p.setData(0.0, 2, 0);
	p.setData(0.0, 3, 0);
}

void Poly3D_dy(Matrix& p, double x, double y, double z) {
	p.setData(0.0, 0, 0);
	p.setData(0.0, 1, 0);
	p.setData(1.0, 2, 0);
	p.setData(0.0, 3, 0);
}

void Poly3D_dz(Matrix& p, double x, double y, double z) {
	p.setData(0.0, 0, 0);
	p.setData(0.0, 1, 0);
	p.setData(0.0, 2, 0);
	p.setData(1.0, 3, 0);
}

//**********************************************************
//forward computation: without epi
//*********************************************************
void FECG(char** filNamBEM, char** filNamMFree, char** filOutput,
		char** TransOutput) {
	matEg engine(MAX_BUF);

	cout << "\n **construct BEM class**\n" << endl;
	BEM_source bem(filNamBEM, filNamMFree, Tri3);

	cout << "\n **BEM Computing**\n" << endl;
	bem.BEMCompute(filOutput);

	int numSur = bem.getNumSur();
	int numVol = bem.getNumVol();
	double* K = bem.getBEMK();
//	mwArray K(numSur,numSur,bem.getBEMK());
//	transpose(K,numSur,);
	double* Ka = bem.getBEMKa();
//	mwArray Ka(1,numSur,bem.getBEMKa());
	double* M = bem.getMFree();
//    transpose(M);

	engine.Open(NULL);
	engine.addData(K, "K", numSur, numSur);
	engine.addData(Ka, "Ka", 1, numSur);
	engine.evalString("K = K.'");
	engine.evalString("K = vertcat(K,Ka)");
	engine.addData(M, "M", numVol, numSur + 1);
	engine.evalString("M = M.'");
//	mwArray M(numVol,numSur+1,bem.getMFree());
//	mwArray Trans = inv(ctranspose(K)*K)*ctranspose(K)*M;
	engine.evalString(" Trans = inv(ctranspose(K)*K)*ctranspose(K)*M");

	/*	FILE* fid = fopen(TransOutput,"wb");
	 double* d = mxGetPr(Trans.GetData());
	 int row,col;
	 row = size(&col,Trans);
	 fwrite(d,sizeof(double),row*col,fid);
	 fclose(fid);
	 */
	cout << "writing Trans files \n" << endl;
	double* Trans = new double[numSur * numVol];
	engine.getData(Trans, "Trans");
	writeFile(TransOutput[0], Trans, sizeof(double), numSur * numVol, "wb");
	delete[] Trans;
	engine.Close();
}

//**********************************************************
//foward computation: with epi
//**********************************************************
void FECG(char** filNamBEM, char** filNamMFree, char** filOutMFree,
		char** filOutBEM, char** TransOutput) {
	matEg* engine = new matEg(MAX_BUF);
	//Laplacian part
	BEM_nosource bem2(filNamBEM, filNamMFree, Tri3);

	bem2.BEMCompute(filOutBEM);
	int numSur_t = bem2.getNumSur();
	double* Kt = bem2.getBEMK();
	double* Pt = bem2.getBEMP();
	for (int i = 0; i < numSur_t * numSur_t; i++)
		Pt[i] = DKe / DTe * Pt[i];

	//Possoin part
	BEM_source bem(filNamBEM, filNamMFree, Tri3);
	bem.BEMCompute(filOutMFree);
	int numSur_h = bem.getNumSur();
	int numVol = bem.getNumVol();
	double* Kh = bem.getBEMK();
	double* Ph = bem.getBEMP();
	double* Ma = bem.getMFree();

	//begin engine***********************
	engine->Open(NULL);
	engine->addData(Kt, "Kt", numSur_t, numSur_t);
	engine->evalString("Kt = transpose(Kt)");
	engine->addData(Pt, "Pt", numSur_t, numSur_t);
	engine->evalString("Pt = transpose(Pt)");
	engine->addData(Kh, "Kh", numSur_h, numSur_h);
	engine->evalString("Kh = transpose(Kh)");
	engine->addData(Ph, "Ph", numSur_h, numSur_h);
	engine->evalString("Ph = transpose(Ph)");
	engine->addData(Ma, "Ma", numVol, numSur_h + 1);
	engine->evalString("numSur_h = size(Kh,1)");
	engine->evalString("numSur_t = size(Kt,1)");
	engine->evalString("numVol = size(Ma,1)");
	engine->evalString("Ma = transpose(Ma)");
	engine->evalString("M = Ma(1:numSur_h,:)");

	//begin to couple
	engine->evalString("K = horzcat(Kh,zeros(numSur_h,numSur_t - numSur_h))");
	engine->evalString("K = vertcat(K,zeros(numSur_t - numSur_h,numSur_t))");

	engine->evalString("P = horzcat(Ph,zeros(numSur_h,numSur_t - numSur_h))");
	engine->evalString(
			"I = horzcat(zeros(numSur_t - numSur_h,numSur_h),eye(numSur_t - numSur_h))");
	engine->evalString("P = vertcat(P,I)");

	engine->evalString("P = P*inv(Pt)*Kt");

	engine->evalString(
			"Trans = inv(K+P)*vertcat(M,zeros(numSur_t - numSur_h,numVol))");

	double* Trans = new double[numSur_t * numVol];
	engine->getData(Trans, "Trans");
	writeFile(TransOutput[0], Trans, sizeof(double), numSur_t * numVol, "wb");
	delete[] Trans;

	engine->Close();
	delete engine;

}

//************************************************
//computeFECG
//***********************************************
int computeFECG() {
	char cont = 'Y';
	do {
		int type;
		cout
				<< "select the forumlation: 1, without heart surface; 2, with heart surface"
				<< endl;
		cin >> type;
		int num_file = 0;
		//ReadMe: filNamBEM**********************************************************
		//with epi: inner suface (cor, seq); outer surface (cor, seq)
		// without epi: surface (cor,seq)
		//***for torso and heart surface, existent in common directory****************
		if (type == 1)
			num_file = 6;
		else
			num_file = 11;
		if (type == 1 || type == 2) {
			char* copath_ref = "/home/wth4280/Documents/Research/data/KIT/Ref/";
			//"Ref/";Auckland/Auckland_f02_2_5/Processed/Results/Ref/";
			//;"../BEMForward/Period_III/Sphere/input/test/";
			char* in = "/home/wth4280/Documents/Research/data/KIT/Final/";
			//"Auckland/Auckland_f02_2_5/Final/";PhysioNet/case2/Auckland/Auckland_f02_2_5/Processed/Results/Final/";
			//"../BEMForward/Period_III/Sphere/input/test/";
			char* out = "/home/wth4280/Documents/Research/data/KIT/Simulation/";
			//"Auckland/Auckland_f02_2_5/Simulation/";Auckland/Auckland_f02_2_5/Processed/Results/Simulation/";
			////"../BEMForward/Period_III/Sphere/output/test/F1/";	PhysioNet/case2
			char id[5];
			sprintf(id, "%d", NUM_NODE);
			char* copath_in = new char[MAX_FILENAME];
			char* copath_out = new char[MAX_FILENAME];
			char** ptr1, **ptr2;
			char** namBEM = new char*[num_file];
			if (type == 1) {
				namBEM[0] = "torso.cor"; //"_tet_vs.cor";//"sphCor_sur_outer_1704.bin"; //"case0003_b352_registered.cor";//"torso_2869.cor"; //inner_1564.bin"; //"torso_new.cor"; //"torso_low.cor";///auckland_sur1254.cor";//".\\cube343\\sphCor.bin";
				namBEM[1] = "torso.tri"; //"peri.tri"; //"sphTri_sur_outer_1704.bin"; //"case0001_b352.seq";//"torso.seq";//inner_1564.bin"; //"torso.seq"; //auckland_sur1254.sur";//".\\cube343\\sphSur.bin";
				namBEM[2] = NULL; //"Forbidden_outer_1704.bin"; //inner_1564.bin";
				namBEM[3] = "torso.center"; //"BEM_element_center_o1704.cor"; //utah.cor";//
				namBEM[4] = "torso_tet.cor";    //"sphCor_sur_outer_1704.bin";//
				namBEM[5] = "torso.tet";    	//"sphTet_sur_outer_1704.bin";//
			} else {
				namBEM[0] = "peri.cor"; //"sphCor_sur_inner_1136.bin";//"auckland_sur1254.cor";//
				namBEM[1] = "peri.tri"; //"sphTri_sur_inner_1136.bin";//"auckland_sur1254.sur";//
				namBEM[2] = NULL; //"Forbidden_inner_1136.bin"; //"torso_new.cor"; //"torso_low.cor";///auckland_sur1254.cor";//".\\cube343\\sphCor.bin";
				namBEM[3] = "peri_tri_center.cor"; //"BEM_element_center_inner_1136.cor";//
				namBEM[4] = "torso.cor"; //"sphCor_sur_outer_68.bin";//"case0003_b352_registered.cor";//"torso_2869.cor";//
				namBEM[5] = "torso.tri"; //"sphTri_sur_outer_68.bin";//"case0001_b352.seq";//"torso.seq";//
				namBEM[6] = NULL; //"Forbidden_outer_68.bin"; //"Forbidden_outer_case0003.bin";
				namBEM[7] = "composite.cor";
				namBEM[8] = "composite.tri";
				namBEM[9] = "Forbidden.bin";
				namBEM[10] = "composite_tri_center"; //"BEM_element_center_whole_case0003.cor";
			}
			char** filNamBEM = new char*[num_file];
			addPath(filNamBEM, in, namBEM, num_file);
			if (num_file == 6)
				addPath(filNamBEM, copath_ref, namBEM, num_file);
			else {
				addPath(filNamBEM, copath_ref, namBEM, 4);
				ptr1 = &filNamBEM[4];
				ptr2 = &namBEM[4];
				addPath(ptr1, in, ptr2, num_file - 4);
			}

			//***for mfree nodes, existent in specific subdirectory********
			strcpy(copath_in, in);
			strcat(copath_in, id);
			strcat(copath_in, "/");
			strcpy(copath_out, out);
			strcat(copath_out, id);
			strcat(copath_out, "/Input/");

			// sourcefile*****************************************************************************
			// the same for both type
			char** filNamMFree = new char*[4];
			char** namMFree = new char*[4];
			namMFree[0] = "heart.cor"; //"sphCor_mfree_1388.bin";//"auckland_1720.cor";//
			namMFree[1] = "heart.cor"; //"sphFib_mfree_1388.bin";// "auckland_1720.fib"; //
			namMFree[2] = "heart_tet.cor"; //"sphCor_sur_inner_1136.bin";// "auckland_sur1254.cor";//"human_sur1911.cor";
			namMFree[3] = "heart.tet"; //"sphTet_sur_inner_1136.bin";//"auckland_sur1254.tet";// "human_sur1911.tet";
			//       	addPath(filNamMFree,in,namMFree,4);
			addPath(filNamMFree, copath_in, namMFree, 2);
			ptr1 = &filNamMFree[2];
			ptr2 = &namMFree[2];
			addPath(ptr1, copath_ref, ptr2, 2);

			//outputfile*******************************************************************************
			//with epi, two group, one is the same with without epi, another is for pure BEM
			char** filOutMFree = new char*[4];
			char** outMFree = new char*[4];
			outMFree[0] = "BEMK.bin";
			outMFree[1] = "BEMP.bin";
			outMFree[2] = "Ka.bin";
			outMFree[3] = "MFree.bin";
			//       addPath(filOutMFree,out,outMFree,4);
			addPath(filOutMFree, copath_out, outMFree, 4);
			char** filOutBEM = new char*[2];
			char** outBEM = new char*[2];
			outBEM[0] = "BEMK_outer.bin";
			outBEM[1] = "BEMP_outer.bin";
//        	addPath(filOutBEM,out,outBEM,2);
			addPath(filOutBEM, copath_out, outBEM, 2);
			//final output for Transfer matrix
			char** finOutput = new char*[1];
			char* output = "Trans.bin";
			//       	addPath(finOutput,out,&output,1);
			addPath(finOutput, copath_out, &output, 1);

			//begin***************************************************************************************
			if (type == 2) {
				cout << "Start FECG" << endl;
				FECG(filNamBEM, filNamMFree, filOutMFree, filOutBEM, finOutput);
				cont = 'N';
			} else {
				FECG(filNamBEM, filNamMFree, filOutMFree, finOutput);
				cont = 'N';
			}

			//***free**********************************
			delete[] copath_in;
			delete[] copath_out;
			for (int i = 0; i < num_file; i++)
				delete[] (filNamBEM[i]);
			delete[] filNamBEM;
			delete[] namBEM;
			for (int i = 0; i < 4; i++)
				delete[] (filNamMFree[i]);
			delete[] filNamMFree;
			delete[] namMFree;
			for (int i = 0; i < 2; i++)
				delete[] (filOutBEM[i]);
			delete[] filOutBEM;
			delete[] outBEM;
			for (int i = 0; i < 4; i++)
				delete[] (filOutMFree[i]);
			delete[] filOutMFree;
			delete[] outMFree;
		} else {
			cout << "no option provided! continue?(y/n)" << endl;
			cin >> cont;
		}
	} while (cont == 'Y' || cont == 'y');

	cout << "ECG Forward Modeling Succeed! ^___^" << endl;
	return EXIT_SUCCESS;
}

//************************************************
//testFECG
//all levels of test, from single element integral, to whole surface BEM, as
//well as MFree method
//************************************************
void testFECG() {
	int type = 0;
	cout << "choose from the following testings: " << endl;
	cout
			<< "1, single elment integral;   2, whole surface integral;   3, MFree integral"
			<< endl;
	cin >> type;
	if (type == 1) {
		BEMCell* bemCell = new Tri3Cell();
		double err = bemCell->test();
		cout << err << endl;
	} else if (type == 2) {
	} else if (type == 3) {
	} else {
		cout << "Warning: no testing functions existent!" << endl;
		exit(0);
	}

}

//***************************************
//check whether a matrix has nan
//****************************************
int checkNan(Matrix& p) {
	double* ptr = p.getDataPtr();
	int num = p.getCol() * p.getRow();
	for (int i = 0; i < num; i++) {
		if (isnan(ptr[i]))
			return 1;
	}
	return 0;
}

