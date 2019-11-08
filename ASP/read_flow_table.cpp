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

struct patch_struct* read_flow_table(struct patch_struct *patch, char *flowtable_in,int *ppatch_num) {

	int itmp;

	//cout << "1.. READING FLOW TABLE:" << endl;

	ifstream FlowTable(flowtable_in, ios::in);

	//3 READ FLOW TABLE
	FlowTable >> *ppatch_num;//skip first patch num value
	int patch_num = *ppatch_num;

	//allco memory
	patch= new struct patch_struct[patch_num]{};

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
			FlowTable >> patch[pch].neigh_ID[neigh];
			FlowTable >> patch[pch].neigh_gamma[neigh];
			patch[pch].pch = pch;
		}


		//init some attributes 
		if (patch[pch].landID == 1)
			patch[pch].channel_acc = 1;//init as 1
		else
			patch[pch].channel_acc = 0;//init as 0

		if(patch[pch].landID==0)
			patch[pch].hill_acc = 1;//init as 1
		else
			patch[pch].hill_acc = 0;//init as 0


		patch[pch].above_num = 0;
		if (patch[pch].landID == 1) 
		patch[pch].longest_stream_length = 1;

	}//

	//compute above information and route channel_acc
	for (int pch = 0; pch != patch_num; pch++) {

		//neighbors
		for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++) {


			for (int pch_next = pch + 1; pch_next != patch_num; pch_next++) {

				if (patch[pch_next].patchID == patch[pch].neigh_ID[neigh]) {//neighbor captured

					//for pch
					patch[pch].neigh_pch[neigh] = pch_next;

					//for neigh
					patch[pch_next].above_pch[patch[pch_next].above_num] = pch;
					patch[pch_next].above_ID[patch[pch_next].above_num] = patch[pch].patchID;
					patch[pch_next].above_gamma[patch[pch_next].above_num] = patch[pch].neigh_gamma[neigh];
					patch[pch_next].above_pch[neigh] = pch;//init as 1

					//be cautious when it excels or equal to 2 means two rivers converges
					patch[pch_next].above_num += 1;//move to the next

					//accumulate channel acc
					patch[pch_next].channel_acc += patch[pch].channel_acc*patch[pch].neigh_gamma[neigh];
					patch[pch_next].hill_acc += patch[pch].hill_acc*patch[pch].neigh_gamma[neigh];

					break;//break out of this for
				}
			}
		}
	}//

	
	for (int pch = 0; pch != patch_num; pch++) {

		//init re_acc_area re_channel_acc
		patch[pch].re_channel_acc = patch[pch].channel_acc;
		patch[pch].re_acc_area = patch[pch].acc_area;

		//setting outlet and converge point
		if (patch[pch].landID == 1) {
			for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++)
				if (patch[(patch[pch].neigh_pch[neigh])].above_num >= 2) {

					patch[pch].outlet_flag = 1;
					patch[(patch[pch].neigh_pch[neigh])].converge_flag = 1;

				}
		}
	}

    //compute stream width fuction
	for (int pch = patch_num-1; pch >= 0; pch--) {
		
		if(patch[pch].landID==1){

		  for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++)
		   patch[pch].longest_stream_length+= patch[(patch[pch].neigh_pch[neigh])].longest_stream_length;

		}
	}
	int total = 0, llength = 0;
	for (int pch = patch_num - 1; pch >= 0; pch--) {
		if (patch[pch].landID == 1) {
			total++;
			if (patch[pch].longest_stream_length > llength) {
				llength = patch[pch].longest_stream_length;
			}
		}
	}
	cout << total / llength << "\n" << llength << "\n" << total << endl;

	//cout << "COMPLETED IN FLOW TABLE\n" << endl;

	return (patch);
}

