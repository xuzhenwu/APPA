#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include "APPA.h"

using namespace std;

struct patch_object* read_flow_table(char* flowtable_in, struct status_information* INF)
{


	cout << "read flow routing information:\n" ;

	ifstream FlowTable(flowtable_in, ios::in);
	if (!FlowTable.is_open()) {
		cout << "unable to open flow routing file: " << flowtable_in << '\n' << endl;
		exit(0);
	}

	FlowTable >> INF->patch_num;//skip first patch num value
	
	int patch_num = INF->patch_num;

	//allco memory
	struct patch_object* patch= new struct patch_object[patch_num]{};
	for (int pch = 0; pch != patch_num; pch++) {

		//read first line parameters
		FlowTable >> patch[pch].patchID;
		FlowTable >> patch[pch].lon;
		FlowTable >> patch[pch].lat;
		FlowTable >> patch[pch].z;
		FlowTable >> patch[pch].acc_area;
		FlowTable >> patch[pch].landID;
		FlowTable >> patch[pch].runoff_coefficient;
		FlowTable >> patch[pch].neigh_num;

		//read neigh ID and gamma
		for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++) {
			FlowTable >> patch[pch].neigh_pch[neigh];
			FlowTable >> patch[pch].neigh_ID[neigh];
			FlowTable >> patch[pch].neigh_gamma[neigh];
			
		}

	}//end of reading flow table
	
	 
	//init attributes 
	for (int pch = 0; pch != patch_num; pch++) {
		
		patch[pch].pch = pch;

		if (patch[pch].landID == 1) {
			patch[pch].channel_acc = 1;//init as 1
			INF->channel_num++;
		}
		else
			patch[pch].channel_acc = 0;//init as 0

		if (patch[pch].landID == 0) {
			patch[pch].hill_acc = 1;//init as 1
			INF->hillslope_num++;
		}
		else
			patch[pch].hill_acc = 0;//init as 0

		//upstream information
		patch[pch].above_num = 0;

		if (patch[pch].landID == 1)
			patch[pch].longest_stream_length = 1;
	
	}//end of init attributes

	//accumulate values via flow routing
	for (int pch = 0; pch != patch_num; pch++) {

		//neighbors
		for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++) {

			int pch_next = patch[pch].neigh_pch[neigh];

		    //accumulated runtime of channel processes
			patch[pch_next].channel_acc += patch[pch].channel_acc*patch[pch].neigh_gamma[neigh];

			//accumulated runtime of hillslope processes
			patch[pch_next].hill_acc += patch[pch].hill_acc*patch[pch].neigh_gamma[neigh];
			
			//upstream information
			patch[pch_next].above_pch[patch[pch_next].above_num] = pch;
			patch[pch_next].above_ID[patch[pch_next].above_num] = patch[pch].patchID;
			patch[pch_next].above_gamma[patch[pch_next].above_num] = patch[pch].neigh_gamma[neigh];
			patch[pch_next].above_pch[neigh] = pch;
			patch[pch_next].above_num += 1;//be cautious when it excels or equal to 2 means two rivers converges

		}
	}//end of accumulation

	//compute converge information
	for (int pch = 0; pch != patch_num; pch++) {

		//setting outlet and converge point
		if (patch[pch].landID == 1) {
			for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++)
				if (patch[(patch[pch].neigh_pch[neigh])].above_num >= 2) {

					patch[pch].outlet_flag = 1;// a pontential outlet partitioned at the outlet 
					patch[(patch[pch].neigh_pch[neigh])].converge_flag = 1;// a potential converge/starting point

				}
		}
	}

    //compute longest_stream_length and estimate theoretical maximum speedup ratio for channel processes
	//accumulate longest_stream_length
	for (int pch = patch_num-1; pch >= 0; pch--) {
		
		if(patch[pch].landID==1){

		  for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++)
		   patch[pch].longest_stream_length+= patch[(patch[pch].neigh_pch[neigh])].longest_stream_length;

		}
	}
	//find maximum length: llength
	double total = 0, llength = 0;
	for (int pch = patch_num - 1; pch >= 0; pch--) {
		if (patch[pch].landID == 1) {
			if (patch[pch].longest_stream_length > llength) {
				llength = patch[pch].longest_stream_length;
			}
			if (patch[pch].channel_acc > total) {
				total = patch[pch].channel_acc;
			}
		}
	}
	//compute theoretical maximum speedup ratio of chanel processes
	INF->TMSRch = total / (llength);

	cout << "complete\n" << endl;
	return (patch);//return the point of the patch array
}

