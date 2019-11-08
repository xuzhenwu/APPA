#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include "PARA.h"

using namespace std;

void print_threads_map(struct patch_struct *patch,char *geo_files,char *threads_out, int patch_num,int thread, int print_flag) {

	//read .patch file 
	//add a tail
	char num[10]{};
	char patch_in[120]{};
	_itoa_s(thread, num, 10);
	char fileout[110]{};

	strcpy_s(fileout, threads_out);
	strcat_s(fileout, num);
	
	strcpy_s(patch_in, geo_files);
	strcat_s(patch_in, ".pch");

	ifstream PatchIn(patch_in, ios::in);
	fstream PrintOut(fileout, ios::out);

	//info
	int  ncols, nrows;
	double xllcorner, yllcorner, cellsize, NODATA_value;
	char header[6][60]{};

	//-----------------------------------------------------------------------------------------------------------------------
	//5. FILES WRINTING 
	//-----------------------------------------------------------------------------------------------------------------------
	//cout << "\n5.. WRINTING FILES\n" << endl;
	//cout <<num<< fileout << endl;
	//read header
	PatchIn >> header[0] >> ncols;
	PatchIn >> header[1] >> nrows;
	PatchIn >> header[2] >> xllcorner;
	PatchIn >> header[3] >> yllcorner;
	PatchIn >> header[4] >> cellsize;
	PatchIn >> header[5] >> NODATA_value;

	//write header
	PrintOut << header[0] << "\t" << ncols << endl;
	PrintOut << header[1] << "\t" << nrows << endl;
	PrintOut << header[2] << "\t" << xllcorner << endl;
	PrintOut << header[3] << "\t" << yllcorner << endl;
	PrintOut << header[4] << "\t" << cellsize << endl;
	PrintOut << header[5] << "\t" << NODATA_value << endl;


	//init new array
	int *pch_ID = new int[ncols*nrows]{};
	int *print_ID = new int[ncols*nrows]{};


	//PRINT THE TABLE
	for (int r = 0; r != nrows; r++) {
		for (int c = 0; c != ncols; c++)
		{
			int inx = r * ncols + c;
			int pchID=0;
			PatchIn >> pch_ID[inx];

			//find patch and it's basinID
			if (pch_ID[inx] >= 0) {
		
					pchID = pch_ID[inx];
					if (print_flag == 1) {
							
						print_ID[inx] = patch[pchID].sthread;
					}
					else{
						print_ID[inx] = patch[pchID].cthread;
						//print_ID[inx] = patch[pchID].channel_acc;
					}
			}
			else print_ID[inx] = -9999;//nodata
			

			//if (inx % 100000 == 0) cout << inx << endl;
			
			//out basinID
			PrintOut << setw(10) << print_ID[inx];
		}

		//change rows
		PrintOut << endl;
	}

	delete[] pch_ID;
	delete[] print_ID;
	PrintOut.close();
	PatchIn.close();

}