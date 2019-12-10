#ifndef  CCreatPitsFile_H_
#define  CCreatPitsFile_H_


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Blender.h"
#include <string>


using namespace std;
void CCreatPitsFile::initialize(){
	for (int i = 0; i < MAXS; i++){
		name1[i] = { NULL };
		name2[i] = { NULL };
	}
};

void CCreatPitsFile::creatPitsFile(char *name1, char*name2, ofstream& out1, ofstream&out2, char *input_prefix)
{
	// open some diagnostic output files

	initialize();

	strcpy(name1, input_prefix);
	strcat(name1, ".build");
	out1.open(name1, ios::out);
	if (!out1){
		cerr << "cannot open build file \n" << endl;
	}

	strcpy(name2, input_prefix);
	strcat(name2, ".pit");
	out2.open(name2, ios::out);
	if (!out2){
		cerr << "cannot open pit file \n" << endl;
		exit(1);
	}

}

#endif