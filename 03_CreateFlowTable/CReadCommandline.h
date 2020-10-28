#ifndef  CReadCommandline_H_
#define  CReadCommandline_H_

#include <iostream>
#include <string>
#include "Blender.h"


void CReadCommandline::initialize(){
	
	CReadCommandline::vflag = 0;		    // verbose flag					 
	CReadCommandline::fl_flag = 0;	    // roads to lowest flna			 
	CReadCommandline::fh_flag = 0;	    // roads to highest flna		 
	CReadCommandline::s_flag = 1;		    // printing stats flag			 
	CReadCommandline::r_flag = 0;		    // road stats flag				 
	CReadCommandline::sc_flag = 1;		// stream connectivity flag		 
	CReadCommandline::st_flag = 0;		// scaling stream side patches 	 
	CReadCommandline::arc_flag = 1;		// arcview input data flag		 
	CReadCommandline::prefix_flag = 0;	// input prefix flag 

	//DEFINE ALGORITHM FLAG
	CReadCommandline::algorithm_flag = 6 ;	//0:D8;1:MD8;2;D_INF;3:RMD_INF;4:MFD-MD;5:MD_inf; 6:parallel/D8
	
	//pch should be true since it's useful for resoting CHESS -p file into images
	CReadCommandline::pch = true;


	//TOPO INDICES
	
	/*
	CReadCommandline::acc = true;//ACC AREA
	CReadCommandline::dts = true;//HORIZONTAL FLOW DISTANCE TO STREAM  
	CReadCommandline::sts = true;//HORIZONTAL FLOW STEPS TO STREAM
	CReadCommandline::dtr = true;//HORIZONTAL FLOW DISTANCE FROM RIDGE 
	CReadCommandline::str = true;//HORIZONTAL FLOW STEPS FROM RIDGE
	CReadCommandline::tga = true;//TOTAL GAMMA (FLOW GENERATION)
	CReadCommandline::twi = true;*/
	
	
	
	CReadCommandline::acc = true;//ACC AREA
	CReadCommandline::dts = true;//HORIZONTAL FLOW DISTANCE TO STREAM  
	CReadCommandline::sts = true;//HORIZONTAL FLOW STEPS TO STREAM
	CReadCommandline::dtr = true;//HORIZONTAL FLOW DISTANCE FROM RIDGE 
	CReadCommandline::str = true;//HORIZONTAL FLOW STEPS FROM RIDGE
	CReadCommandline::tga = true;//TOTAL GAMMA (FLOW GENERATION)
	CReadCommandline::twi = true;



	if (CReadCommandline::algorithm_flag == 6) {
		CReadCommandline::acc = false;
		CReadCommandline::dts = false;
		CReadCommandline::sts = false;
		CReadCommandline::dtr = false;
		CReadCommandline::str = false;
		CReadCommandline::tga = false;
		CReadCommandline::twi = false;
	}

	
	CReadCommandline::width = 5;			// default road width
	input_prefix = new char [140];
}

void CReadCommandline::Commandline(int argc, char *argv[])
{
	initialize();

	int i = 1;
	while (i < argc) {

		if (strcmp(argv[i], "-sc") == 0){
			i += 1;
			if (i == argc)  {
				std::cerr << "\n FATAL ERROR: Output prefix not specified \n" << std::endl;
				exit(0);
			}
			sc_flag = (int)atoi(argv[i]);
		}

		if (strcmp(argv[i], "-s") == 0){
			s_flag = 1;
		}

		if (strcmp(argv[i], "-w") == 0){
			i += 1;
			if (i == argc)  {
				std::cerr << "\n FATAL ERROR: Output prefix not specified \n" << std::endl;
				exit(0);
			}
			width = (double)atof(argv[i]);
		}

		if (strcmp(argv[i], "-res") == 0){
			i += 1;
			if (i == argc)  {
				std::cerr << "\n FATAL ERROR: Data resolution not specified \n" << std::endl;
				exit(0);
			}
			cell = (double)atof(argv[i]);
		}

		if (strcmp(argv[i], "-pre") == 0){
			i += 1;
			if (i == argc)  {
				std::cerr << "\n FATAL ERROR: Output prefix not specified\n" << std::endl;
				exit(0);
			}

			//-------------------------------------------------------------------------//
			// Allocate an array for the output prefix and read in the output prefix   //
			//-------------------------------------------------------------------------//
			input_prefix = new char[1 + strlen(argv[i])];
			if (!input_prefix){
				std::cerr << "\n FATAL ERROR: Cannot allocat output_prefix\n" << std::endl;
				exit(0);
			}
			strcpy(input_prefix, argv[i]);
			prefix_flag = 1;
		}

		if (strcmp(argv[i], "-r") == 0){
			r_flag = 1;
		}

		if (strcmp(argv[i], "-a") == 0){
			arc_flag = 1;
		}
		i += 1;
	} //end of while

	//if (prefix_flag == 0) {
		//std::cerr << "\n FATAL ERROR: You must specify a prefix for image files\n" << std::endl;
		//exit(0);
	//}

	if ((fl_flag) || (fh_flag))
		f_flag = 1;
	else
		f_flag = 0;

}

#endif