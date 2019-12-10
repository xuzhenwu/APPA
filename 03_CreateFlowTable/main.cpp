//==============================================================================================
//  Program name: create_flowpaths								
//  Options                                                     
//                                                              
//		-v 	Verbose Option										
//		-fl roads to lowest flna interval 					    
//		-fh roads to highest flna interval					    
//		-s  print drainage statistics							
//		-r 	road flag for drainage statistics					
//		-sc	stream connectivity is assumed						
//			1	random slope value								
//			2	internal slope value							
//			0	no connectivity (default)						
//		-sd scale dem values by this amount
//		-st	scale streamside transmissivity						
//		-pre	input image file name prefix					
//		-w 	road width (default is 5m)							
//		-a arcview ascii data files (default is GRASS ascii)	
//                                                              
//  DESCRIPTION                                                 
//=============================================================================================

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <thread>
#include "Blender.h"
#include "CReadArcInfoImage.h"
#include "CReadCommandline.h"
#include "CreatePitsFile.h"
#include "CreatFlowTable.h"

using std::thread;

int main(int argc, char *argv[])
{

	CReadCommandline *pCommandline = new CReadCommandline;
	CReadArcInfoImages *pCArcImages = new CReadArcInfoImages;
	CCreatPitsFile  *pPitsFile = new CCreatPitsFile;

	//Read commandline arguments
	pCommandline->Commandline(argc, argv);
	
	//initialize paramters for ArcInfo images
	pCArcImages->initialize();

	//Assign ArcInfo image names amd read images
	pCArcImages->input_prompt(pCArcImages->maxr, pCArcImages->maxc, pCommandline->input_prefix,
		pCArcImages->fndem, pCArcImages->fnslope, pCArcImages->fnpatch, pCArcImages->fnstream,
		pCArcImages->fnroads, pCArcImages->fntable,pCArcImages->fnsubbasin, pCommandline->arc_flag);

	//Creat two files to record pits
	pPitsFile->creatPitsFile(pPitsFile->name1, pPitsFile->name2, pPitsFile->out1, pPitsFile->out2, pCommandline->input_prefix);
	
	//allocate flow table
	CreatFlowTable *pCreatFlowTable = new class CreatFlowTable;
	/*pCreatFlowTable->initialize();*/
	flow_struct *flow_table_ = new flow_struct[(pCArcImages->maxr)*(pCArcImages->maxc)]{};
	pCreatFlowTable->flow_table = flow_table_;

	//sorting the flow table in an desending order
	pCreatFlowTable->num_patches = pCreatFlowTable->sort_flow_table(pCArcImages->patch_pch, pCArcImages->dem, pCArcImages->patch, pCArcImages->maxr, pCArcImages->maxc );

	//init all variables
	pCreatFlowTable->zero_flow_table(pCreatFlowTable->flow_table, pCArcImages->maxr, pCArcImages->maxc);

	
	//regouping for parallel 
	for (int thread_inx=0; thread_inx != THREADS +1; thread_inx++) {
		
		pCommandline->rcboun[thread_inx] = int((pCArcImages->maxr * pCArcImages->maxc /(THREADS * 1.0))*(thread_inx));
		pCommandline->pchboun[thread_inx] =(int) (( pCreatFlowTable->num_patches/ (THREADS*1.0))*(thread_inx));

	}

	//build the structure of flow tables
	thread build_thd[THREADS];
	for (int thread_inx=0; thread_inx != THREADS; thread_inx++){
		
		build_thd[thread_inx]= thread(build_flow_table,pCreatFlowTable->flow_table, pCArcImages->dem, pCArcImages->slope, pCArcImages->patch, pCArcImages->stream, pCArcImages->roads, pCArcImages->subbasin, pCArcImages->flna, pCArcImages->maxr,
		pCArcImages->maxc, pCommandline->f_flag, pCommandline->sc_flag, pCArcImages->cell, pCArcImages->xllcorner, pCArcImages->yllcorner,
			pCArcImages->patch_pch,pCommandline->rcboun,thread_inx);
	}
	for (int thread_inx=0; thread_inx !=THREADS; thread_inx++) {
		build_thd[thread_inx].join();
	}	

	pPitsFile->out1.close();

	//processes patches - computing means and neighbour slopes and gammas
	thread gamma_thd[THREADS];
	for (int thread_inx = 0; thread_inx != THREADS; thread_inx++) {
		gamma_thd[thread_inx] = thread(compute_gamma,pCreatFlowTable->flow_table, pCreatFlowTable->num_patches, pCommandline->sc_flag,
			pCArcImages->cell, pCommandline->algorithm_flag,pCommandline->pchboun, thread_inx);
	}
	for (int thread_inx = 0; thread_inx != THREADS; thread_inx++) {
		gamma_thd[thread_inx].join();
	}


	// remove pits and re-order patches appropriately
	pCreatFlowTable->remove_pits(pCreatFlowTable->flow_table, pCreatFlowTable->num_patches, pCommandline->sc_flag, pCArcImages->cell, pPitsFile->out2);

	// add roads
	pCreatFlowTable->add_roads(pCreatFlowTable->flow_table, pCreatFlowTable->num_patches, pCArcImages->cell);

	// find_receiving patch for flna options
	if (pCommandline->f_flag)
		pCreatFlowTable->route_roads_to_patches(pCreatFlowTable->flow_table, pCreatFlowTable->num_patches, pCArcImages->cell, pCommandline->fl_flag);

	if (pCommandline->s_flag == 1) {
		pCreatFlowTable->compute_upslope_area(pCreatFlowTable->flow_table, pCreatFlowTable->num_patches, pCommandline->r_flag, pCArcImages->cell);
		//print_drain_stats(num_patches, flow_table);
		//compute_dist_from_road(flow_table, num_patches, out2, cell);	
		//compute_drainage_density(flow_table, num_patches, cell);
	}

	
	pCreatFlowTable->print_flow_table(pCreatFlowTable->num_patches, pCreatFlowTable->flow_table, pCommandline->sc_flag, pCArcImages->cell, pCommandline->input_prefix, pCommandline->width, pCommandline->algorithm_flag);

	pCreatFlowTable->print_topo_indices(pCreatFlowTable->num_patches, pCreatFlowTable->flow_table, pCArcImages->patch, pCommandline, pCArcImages->cell, pCommandline->input_prefix, pCommandline->width,
		pCArcImages->maxr, pCArcImages->maxc, pCommandline->algorithm_flag);

	//close output files
	pPitsFile->out2.close();

	//release dynamically allocatable array
	delete[] pCArcImages->patch;
	delete[] pCArcImages->dem;
	delete[] pCArcImages->slope;
	delete[] pCArcImages->stream;
	delete[] pCArcImages->roads;
	if (pCArcImages->flna)
		delete[] pCArcImages->flna;
	delete[] pCreatFlowTable->flow_table;

	cout << "END OF FLOW ROUTING COMPUTATIONS:: CHECK FILES" << endl;
	getchar();
	return 0;
}