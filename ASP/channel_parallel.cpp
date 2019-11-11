//---------------------------------------------------------------------------------------------------------------------------
//MAJOR FUNCTION OF ASP
//---------------------------------------------------------------------------------------------------------------------------
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

void channel_parallel(struct patch_struct *patch, double deviation_index, int threads, struct inf* INF,int out_flag) {

	//=======================================================================================================================
	//1. LOCAL VARS(JUST SKIP IF YOU DONT WANT TO rest_numAD ABOUT IT)
	//=======================================================================================================================

	const int infinite = 99999999;
	// while partitioning, num of subchannels will be more than the channel_num we want
	// thus we need more memory for storage and by using optimize method we finally got sufficient channel_num

	//=======================================================================================================================
	//1. STEP 1:ENIVIRONMENT SETTING
	//=======================================================================================================================
	//1.1 compute INF information
	double smallest_unit = infinite;
	INF->channel_num = 0;
	INF->partitioned_num = 0;
	INF->rest_num = 0;
	for (int pch = 0; pch != INF->patch_num; pch++) {

		//the backup are used for renew
		if (patch[pch].landID == 1) {
			INF->channel_num++;
			patch[pch].re_channel_acc_backup = patch[pch].re_channel_acc;

			if (patch[pch].re_channel_acc > 0 && patch[pch].channel_inx == 0)
				if (patch[pch].outlet_flag == 1) {
					if (patch[pch].re_channel_acc < smallest_unit)
						smallest_unit = patch[pch].re_channel_acc;
				}
			if (patch[pch].channel_state != 0)
				INF->partitioned_num++;
			else
				INF->rest_num++;
		}
	}

	//1.2 allocate memory for threads 
	int **thread_pches = new int *[ threads];
	for (int i = 0; i <  threads; i++)
	{
		thread_pches[i] = new int[INF->patch_num]{};
	}
	int *thread_pch_num = new int[threads] {};

	//1.3 specify objectives: TAI,TAM,TCHM
	//add this line for avoiding re_channel_acc==1 
	smallest_unit = 1;
	double Tchm = smallest_unit * (1 + deviation_index);// Target ACC Arest_numA
	if (Tchm < 1) Tchm = 1;

	//rest of channel grids in a thread
	for (int inx = 0; inx != threads; inx++) {
		INF->TAi[inx] = Tchm;
	}
	double TAi_smallest;
	double TAM;

	INF->channel_unit_num = 0;

	// the ID of channel num 
	int thread_inx = 0;//start from 0 and controls the num of subchannels (before we optimize it)

	//=======================================================================================================================
	//2. STEP 2: SCHEME GENERATION
	//=======================================================================================================================
	int partition_flag;
	//START PARTITIONING
	do {
		//init local parameters in each round of partition
		partition_flag = 0;
		int best_channel_pch_num{};
		int best_channel_pch_ID{};
		double best_deviation=infinite;

		//---------------------------------------------------------------------------------------------------------------------------
		// 2.1 CHANNEL PARTITION
		//---------------------------------------------------------------------------------------------------------------------------

		//define TAM and thread_inx
		TAM = 0;
		for (int inx = 0; inx != threads; inx++) {
			if (INF->TAi[inx] > TAM)
			{
				TAM = INF->TAi[inx]; thread_inx = inx;
			}
		}
		if (fabs(TAM) < 0.000001)//when there is no avaible vacumn
			break;

		//search for satisfied outlet
		for (int pch = 0; pch != INF->patch_num; pch++) {

			//stream patches that are remained 
			if(patch[pch].landID==1)
			if (patch[pch].re_channel_acc > 0)
				if (patch[pch].channel_state == 0)//unmatched 
					if (patch[pch].lock_state == 0)//unlocked downstream points
					{
						//find best ID
						if (fabs(patch[pch].re_channel_acc - TAM) < best_deviation && patch[pch].re_channel_acc <= TAM) {//when there is good partition

							best_channel_pch_num = thread_pch_num[thread_inx];
							best_channel_pch_ID = pch;
							best_deviation = fabs(patch[pch].re_channel_acc - TAM);
							partition_flag = 1;
						}
					}
		}//end of search
		if (partition_flag == 0)// with no satisfactory outlet
			break; //break out of partition
	
		//---------------------------------------------------------------------------------------------------------------------------
		//2.2 THREAD ALLOCATION
		//---------------------------------------------------------------------------------------------------------------------------
		int partition_num = 0;


		//add first ID in the list
		thread_pches[thread_inx][(thread_pch_num[thread_inx])] = patch[best_channel_pch_ID].patchID;

		//add boundary
		thread_pch_num[thread_inx]++;
		partition_num++;

		//define upstream grids
		for (int up_pch_inx = best_channel_pch_ID - 1; up_pch_inx >= 0; up_pch_inx--) {

			// only when it's unlabelled
			if (patch[up_pch_inx].channel_inx == 0 && patch[up_pch_inx].landID == 1) {
				//check it's neigh
				for (int neigh = 0; neigh != patch[up_pch_inx].neigh_num; neigh++) {

					//check the list now in this channels
					for (int channel_pch_inx = 0; channel_pch_inx != thread_pch_num[thread_inx]; channel_pch_inx++) {

						//when it's neigh is in the channel list
						if (patch[up_pch_inx].neigh_ID[neigh] == thread_pches[thread_inx][channel_pch_inx])
							if (patch[up_pch_inx].channel_state == 0)
							{

								//add label
								patch[up_pch_inx].channel_inx = thread_inx + 1;//start from 1
								patch[up_pch_inx].channel_state = 1;//transient state
								patch[up_pch_inx].re_channel_acc = 0;

								//add it's ID in the channel list
								thread_pches[thread_inx][(thread_pch_num[thread_inx])] = patch[up_pch_inx].patchID;
								thread_pch_num[thread_inx]++;
								partition_num++;

								break;
							}
					}
				}
			}
		}//end of seaching it's upslope area


		if (fabs(partition_num - patch[best_channel_pch_ID].re_channel_acc) >0.001)
			cout << "it the channel grids with no eqaul tg or there is an error in partition?" << endl;

		patch[best_channel_pch_ID].channel_inx = thread_inx + 1;//start from 1
		patch[best_channel_pch_ID].channel_state = 1;
		INF->TAi[thread_inx] -= patch[best_channel_pch_ID].re_channel_acc;

		//renew re_channel_pch
		// in an desending order to to final channel outlet
		int up_pch = 0;
		int down_pch = best_channel_pch_ID;
		//searching downslope area
		while (down_pch != INF->patch_num - 1) {

			up_pch = down_pch;

			//search for down pch

			for (down_pch = up_pch + 1; down_pch != INF->patch_num; down_pch++) {
				if (patch[down_pch].patchID == patch[up_pch].neigh_ID[0])
					break;
			}

			patch[down_pch].re_channel_acc -= patch[best_channel_pch_ID].re_channel_acc;
			patch[down_pch].lock_state = 1;
		} //till outlet

		patch[best_channel_pch_ID].re_channel_acc = 0;
		
		INF->channel_unit_num += 1;//an extra unit partitioned

	} while (partition_flag == 1);

	//correlating final channels with final thread
	int *channel_num = new int[threads] {};
	for (int pch = 0; pch != INF->patch_num; pch++) {
		if (patch[pch].channel_state == 1) {
			//start from 1
			patch[pch].layer_thread = patch[pch].channel_inx;
			channel_num[patch[pch].channel_inx - 1]++;
		}
	}

	//compute hsr
	TAM = 0;
	double TV = 0;

	//search for max
	for (thread_inx = 0; thread_inx != threads; thread_inx++) {
		INF->TAi[thread_inx] = Tchm - INF->TAi[thread_inx];
		TV += INF->TAi[thread_inx];
		if (INF->TAi[thread_inx] > TAM)
			TAM = INF->TAi[thread_inx];
	}
	INF->SR=INF->CSR = TV /TAM;

	//might be deleted in time of next check
	INF->rest_num = 0;
	for (int pch = 0; pch != INF->patch_num; pch++) {
		if (patch[pch].landID == 1 && patch[pch].channel_state == 0)
			INF->rest_num++;
	}
	INF->partitioned_num = INF->channel_num - INF->rest_num - INF->partitioned_num;

	//free memory
	for (int i = 0; i < threads; ++i)
		delete[] thread_pches[i];
	delete[] thread_pches;
	delete[] thread_pch_num;

	return;
}//end of channel parallelization