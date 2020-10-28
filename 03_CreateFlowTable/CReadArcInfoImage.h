#ifndef  CReadArcInfoImage_H_
#define  CReadArcInfoImage_H_

#include <iostream>
#include <fstream>
#include "blender.h"


using namespace std;

void CReadArcInfoImages::initialize(){
	fntable[MAXS] = { NULL };
	fnpatch[MAXS] = { NULL };
	fnslope[MAXS] = { NULL };
	fnstream[MAXS] = { NULL };
	fnroads[MAXS] = { NULL };
	fnflna[MAXS] = { NULL };
	//dem = { nullptr }, slope = { nullptr }, flna = { nullptr };
	//patch = { nullptr }, stream = { nullptr }, roads = { nullptr };
}
//--------------------------------------------------------------------------
//	input_header() - input information (row, col) from [root].header
//--------------------------------------------------------------------------
void CReadArcInfoImages::input_header(int *maxr, int *maxc, char *fndem, int arc_flag)
{
	ifstream inFileDem;
	char ChCols[20], ChRows[20];
	char Xllcorner[20], Yllcorner[20];
	char Cell[20];
	double xllcorner, yllcorner;
	double cell;

	inFileDem.open(fndem, ios::in);
	if (inFileDem.is_open()){
		inFileDem >> ChCols >> *maxc;
		//std::cout << " cols are " << *maxc << std::endl;
		CReadArcInfoImages::maxc = *maxc;

		inFileDem >> ChRows >> *maxr;
		//std::cout << " rows are " << *maxr << std::endl;
		CReadArcInfoImages::maxr = *maxr;

		inFileDem >> Xllcorner >> xllcorner;
		CReadArcInfoImages::xllcorner = xllcorner;
		inFileDem >> Yllcorner >> yllcorner;
		CReadArcInfoImages::yllcorner = yllcorner;

		inFileDem >> Cell >> cell;
		CReadArcInfoImages::cell = cell;

		inFileDem.close();
	}
	else {
		std::cout << "cannot open ArcInfo DEM files \n" << std::endl;
	}
}

//-------------------------------------------------------------------------------------
//	input_prompt() - input root filename, create full filenames
//-------------------------------------------------------------------------------------
void	CReadArcInfoImages::input_prompt(int maxr, int maxc, char *filename, char *fndem, char *fnslope, char *fnpatch, char *fnstream, char *fnroads, char *fntable,char *fnsubbasin, int arc_flag)
{
	// copy the root filename into the specific filenames
	strcpy(fndem, filename);
	strcpy(fnslope, filename);
	strcpy(fnpatch, filename);
	strcpy(fnstream, filename);
	strcpy(fnroads, filename);
	strcpy(fntable, filename);
	strcpy(fnsubbasin, filename);

	cout << "read .patch, .dem, .stream file:";

	// append '.' extensions to each filename (these should be generalized)
	strcat(fndem, ".dem");
	strcat(fnslope, ".slope");
	strcat(fnpatch, ".patch");
	strcat(fnsubbasin, ".subbasin");


	//xu. changed in this part
	strcat(fnstream, ".stream");
	strcat(fnroads, ".road");
	strcat(fntable, ".table");

	input_header(&maxr, &maxc, fndem, arc_flag);

	// Dynamically allocate memory for input map images
	patch = new int[maxr*maxc]{};
	input_ascii_int(patch, fnpatch, maxr, maxc, arc_flag);
	CReadArcInfoImages::patch = patch;

	dem = new double[maxr*maxc]{};
	input_ascii_double(dem, fndem, maxr, maxc, arc_flag);
	CReadArcInfoImages::dem = dem;

	slope = new double[maxr*maxc]{};
	//input_ascii_double(slope, fnslope, maxr, maxc, arc_flag);
	CReadArcInfoImages::slope = slope;

	stream = new int[maxr*maxc]{};
	input_ascii_int(stream, fnstream, maxr, maxc, arc_flag);
	CReadArcInfoImages::stream = stream;

	roads = new int[maxr*maxc]{};
	//input_ascii_int(roads, fnroads, maxr, maxc, arc_flag);
	CReadArcInfoImages::roads = roads;

	subbasin = new int[maxr*maxc]{};
	//input_ascii_int(subbasin, fnsubbasin, maxr, maxc, arc_flag);
	CReadArcInfoImages::subbasin = subbasin;


	//init as zero and then it will computed in spr
	CReadArcInfoImages::patch_pch = new int[maxr*maxc]{};//defines the flow_table list order in 1:maxr*maxc

	//if (f_flag) {
	//	flna = new double[*maxr**maxc];
	//	input_ascii_double(flna, fnflna, *maxr, *maxc, arc_flag);
	//}
	//else flna = NULL;
	cout << "complete\n" << endl;
	return;
}

//-------------------------------------------------------------------------
//input_ascii_int() - input an ascii image into an interger array using the
//      (row, col) coordinates maxr and maxc.
//-------------------------------------------------------------------------
void CReadArcInfoImages::input_ascii_int(int *aray, char *filename, int mc, int mr, int arc_flag)
{
	ifstream InPutImage;
	int  r;
	int max;
	char STR[70];
	max = 0;

	InPutImage.open(filename, ios::in);
	if (InPutImage.is_open()){
		// skip header
		if (arc_flag == 0)
		for (r = 0; r < LEN_GRASSHEADER; r++)  InPutImage >> STR;
		else
		for (r = 0; r < LEN_ARCHEADER; r++)  InPutImage >> STR;

		for (r = 0; r < mr*mc; r++) {
			InPutImage >> aray[r];
			if (aray[r] > max)  max = aray[r];
		}

		//std::cout << "\n Max for " << filename << " is " << max << std::endl;
		InPutImage.close();
	}
	else{
		std::cout << "cannot open file " << filename << std::endl;
		exit(0);
	}

	return;
}


//---------------------------------------------------------------------------------------------
// input_ascii_double() - input an ascii image into an double array using the
//                        (row, col) coordinates maxr and maxc.
//---------------------------------------------------------------------------------------------
void CReadArcInfoImages::input_ascii_double(double *aray, char *filename, int mc, int mr, int arc_flag)
{
	ifstream InPutImage;
	int  r;
	double max;
	char STR[70];

	max = 0;
	InPutImage.open(filename, ios::in);
	if (InPutImage.is_open()){
		/* skip header */
		if (arc_flag == 0)
		for (r = 0; r < LEN_GRASSHEADER; r++)  InPutImage >> STR;
		else
		for (r = 0; r < LEN_ARCHEADER; r++)  InPutImage >> STR;

		for (r = 0; r < mr*mc; r++) {
			InPutImage >> aray[r];
			aray[r] = (double)(aray[r]);
			if (aray[r] > max)  max = aray[r];
		}

		//std::cout << "\n Max for " << filename << " is " << max << std::endl;
		InPutImage.close();
	}
	else {
		std::cout << "cannot open file " << filename << std::endl;
		exit(0);
	}

	return;
}

#endif
//---------------------------------------------------------------------------------------------