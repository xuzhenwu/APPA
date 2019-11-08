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

void channel_parallel(struct patch_struct *patch, double Target_Deviation, int Division_Times, int threads, int patch_num, double Resolution, double smallest_portion_channel_unit, struct inf* Out_Inf,int out_flag) {

	//=======================================================================================================================
	//1. LOCAL VARS(JUST SKIP IF YOU DONT WANT TO READ ABOUT IT)
	//=======================================================================================================================

	const int infinite = 99999999;
	// while partitioning, num of subchannels will be more than the channel_num we want
	// thus we need more memory for storage and by using optimize method we finally got sufficient channel_num
	const int sub_channel_times(400);

	int itmp{};
	int precision{};
	int first_round{};
	int round_flag{};
	int round_total{};
	int optimize_flag{};
	int new_main_channel_flag{};
	int best_subchannel_pch_num{};
	int best_subchannel_pch_ID{};
	double best_deviation{};
	int min_subchannel_pch_num{};
	int min_subchannel_pch_ID{};
	int max_flag{};
	int max_channel_ID{};
	int max_subchannel_ID_OUT{};
	int max_channel_pch_num{};
	int optimize_total{};
	int break_flag{};

	//double stream_boun = Stream_Boun_Area * 1000 * 1000 / (Resolution*Resolution);




	//SUB-BASIN_ID IN FINAL THREADS
	int *channel_inx_ID = new int[sub_channel_times * threads]{};
	int *channel_pch_num = new int[sub_channel_times * threads]{};

	//ID of subchannels in a channel 
	int **final_subchannel_collection = new int *[threads] {};
	for (int i = 0; i < threads; ++i)
	{
		final_subchannel_collection[i] = new int[threads*sub_channel_times]{};
	}
	//num of patches in a subchannel 
	int **final_subchannel_pch_num = new int *[threads] {};
	for (int i = 0; i < threads; ++i)
	{
		final_subchannel_pch_num[i] = new int[threads*sub_channel_times]{};
	}

	int *final_channel_pch_num = new int [threads] {};
	//num of patches in a channel
	int *final_pch_num = new int [threads] {};
	//num of subchannels in a channel
	int *final_subchannel_num = new int [threads] {};


	//=======================================================================================================================
	//1. STEP 1:ENIVIRONMENT SETTING
	//=======================================================================================================================
	double smallest_unit = infinite;
	Out_Inf->CN = 0;
	Out_Inf->PN = 0;
	Out_Inf->RE = 0;
	for (int pch = 0; pch != patch_num; pch++) {

		//the backup are used for renew
		if (patch[pch].landID == 1) {

			Out_Inf->CN++;

			patch[pch].re_channel_acc_backup = patch[pch].re_channel_acc;

			if (patch[pch].re_channel_acc > 0 && patch[pch].channel_inx == 0)
				if (patch[pch].outlet_flag == 1) {
					if (patch[pch].re_channel_acc < smallest_unit)
						smallest_unit = patch[pch].re_channel_acc;
				}
			if (patch[pch].channel_state != 0)
				Out_Inf->PN++;
			else
				Out_Inf->RE++;
		}
	}


	//SUB_BASINS CAUGHT IN STEP 1 
	int **channel_pches = new int *[sub_channel_times * threads];
	for (int i = 0; i < sub_channel_times * threads; i++)
	{
		channel_pches[i] = new int[Out_Inf->CN]{};
	}


	//add this line for avoiding re_channel_acc==1 
	smallest_unit = Out_Inf->RE / threads * smallest_portion_channel_unit;
	smallest_unit = 1;
	double AT = smallest_unit * (1 + Target_Deviation);// Target ACC AREA
	if (AT < 1) AT = 1;
	if (Target_Deviation > 20)
		int a = 0;

	//rest of channel grids in a thread
	for (int inx = 0; inx != threads; inx++) {
		Out_Inf->RE_CH[inx] = AT;
	}
	double RE_CH_smallest;
	Out_Inf->CUN = 0;

	// the ID of channel num 
	int channel_inx = 0;//start from 0 and controls the num of subchannels (before we optimize it)

	//=======================================================================================================================
	//2. STEP 2:CHANNEL UNITS PARTITION
	//=======================================================================================================================
	//rest patches waited for distribution
	int partition_flag;
	//START PARTITIONING
	do {

		partition_flag = 0;

		new_main_channel_flag = 0;

		//INIT
		best_subchannel_pch_num = 0;
		best_subchannel_pch_ID = 0;
		best_deviation = infinite;

		//---------------------------------------------------------------------------------------------------------------------------
		// START SEACHING FOR IDEAL SUBBAISN ASSEMBLE FROM THE HIGHEST PATCH 
		//---------------------------------------------------------------------------------------------------------------------------
		for (int pch = 0; pch != patch_num; pch++) {

			//stream patches that are remained 
			if (patch[pch].re_channel_acc > 0)
				if (patch[pch].channel_state == 0)//unmatched 
					if (patch[pch].lock_state == 0)//unlocked downstream points
					{
						//find best ID
						//there is no need for && patch[pch].re_channel_acc <= AT
						if (fabs(patch[pch].re_channel_acc - AT) < best_deviation && patch[pch].re_channel_acc <= AT) {//when there is good partition

							best_subchannel_pch_num = channel_pch_num[channel_inx];
							best_subchannel_pch_ID = pch;
							best_deviation = fabs(patch[pch].re_channel_acc - AT);
							partition_flag = 1;
						}
					}

		}
		if (partition_flag == 0) { // when no satisfactory outlet remained

			break;  //break out of partition
		}
		//END OF SERCHING  

		//---------------------------------------------------------------------------------------------------------------------------
		//STARTING DISTRIBUTION AND RENEW re_channel_acc
		//---------------------------------------------------------------------------------------------------------------------------

		int partition_num = 0;


		//add first ID in the list
		channel_pches[channel_inx][(channel_pch_num[channel_inx])] = patch[best_subchannel_pch_ID].patchID;

		//add boundary
		channel_pch_num[channel_inx]++;
		partition_num++;

		//searching for upslope undefined patches
		for (int up_pch_inx = best_subchannel_pch_ID - 1; up_pch_inx >= 0; up_pch_inx--) {

			// only when it's unlabelled
			if (patch[up_pch_inx].channel_inx == 0 && patch[up_pch_inx].landID == 1) {
				//check it's neigh
				for (int neigh = 0; neigh != patch[up_pch_inx].neigh_num; neigh++) {

					//check the list now in this channels
					for (int channel_pch_inx = 0; channel_pch_inx != channel_pch_num[channel_inx]; channel_pch_inx++) {

						//when it's neigh is in the channel list
						if (patch[up_pch_inx].neigh_ID[neigh] == channel_pches[channel_inx][channel_pch_inx])
							if (patch[up_pch_inx].channel_state == 0)
							{

								//add label
								patch[up_pch_inx].channel_inx = channel_inx + 1;//start from 1
								patch[up_pch_inx].channel_state = 1;//transient state
								patch[up_pch_inx].re_channel_acc = 0;

								//add it's ID in the channel list
								channel_pches[channel_inx][(channel_pch_num[channel_inx])] = patch[up_pch_inx].patchID;
								channel_pch_num[channel_inx]++;
								partition_num++;
								Out_Inf->RE_CH[channel_inx] --;

								break;
							}
					}

				}
			}

		}//end of seaching it's upslope area


		if (partition_num != patch[best_subchannel_pch_ID].re_channel_acc)
			exit(-55);

		patch[best_subchannel_pch_ID].channel_inx = channel_inx + 1;//start from 1
		patch[best_subchannel_pch_ID].channel_state = 1;
		Out_Inf->RE_CH[channel_inx] --;



		//renew re_channel_pch
		// in an desending order to to final channel outlet
		int up_pch = 0;
		int down_pch = best_subchannel_pch_ID;
		//searching downslope area
		while (down_pch != patch_num - 1) {

			up_pch = down_pch;

			//search for down pch

			for (down_pch = up_pch + 1; down_pch != patch_num; down_pch++) {
				if (patch[down_pch].patchID == patch[up_pch].neigh_ID[0])
					break;
			}

			patch[down_pch].re_channel_acc -= patch[best_subchannel_pch_ID].re_channel_acc;
			patch[down_pch].lock_state = 1;
		} //till outlet

		patch[best_subchannel_pch_ID].re_channel_acc = 0;
		
		Out_Inf->CUN += 1;

		//after divided
		channel_inx++;
		//cout << channel_inx << "\t" << channel_pch_num[channel_inx - 1] << "\t" << rest_num << endl;
		//---------------------------------------------------------------------------------------------------------------------------
		// END OF COMPASATION PROCESS FOR FINDING SUBBASINS
		//---------------------------------------------------------------------------------------------------------------------------

		//when there is too much channel units
		RE_CH_smallest = 0;
		for (int inx = 0; inx != threads; inx++) {
			if (Out_Inf->RE_CH[inx] > RE_CH_smallest)
			{
				RE_CH_smallest = Out_Inf->RE_CH[inx];
				channel_inx = inx;
			}
		}
		AT = Out_Inf->RE_CH[channel_inx];


	} while (partition_flag == 1);


	

	//correlating final channels with final thread
	int *channel_num = new int[threads] {};
	for (int pch = 0; pch != patch_num; pch++) {

		if (patch[pch].channel_state == 1) {

			//start from 1
			patch[pch].layer_thread = patch[pch].channel_inx;
			channel_num[patch[pch].channel_inx - 1]++;
		}
	}

	//CHECK SPEED UP RATIO
	double Max_Channel_Thread = 0;
	double Final_Channel_Num = 0;
	//sum up 
	for (channel_inx = 0; channel_inx != threads; channel_inx++) {
		Final_Channel_Num += channel_num[channel_inx];
	}
	//search for max
	for (channel_inx = 0; channel_inx != threads; channel_inx++) {
		if (channel_num[channel_inx] > Max_Channel_Thread)
			Max_Channel_Thread = channel_num[channel_inx];
	}
	Out_Inf->SR = Final_Channel_Num / Max_Channel_Thread;

	Out_Inf->RE = 0;
	for (int pch = 0; pch != patch_num; pch++) {
		if (patch[pch].landID == 1 && patch[pch].channel_state == 0)
			Out_Inf->RE++;
	}
	Out_Inf->PN = Out_Inf->CN - Out_Inf->RE - Out_Inf->PN;


	//FREE MEMORY
	for (int i = 0; i < sub_channel_times * threads; ++i)
		delete[] channel_pches[i];
	delete[] channel_pches;

	delete[] channel_inx_ID;
	delete[] channel_pch_num;

	for (int i = 0; i < threads; ++i)
		delete[] final_subchannel_collection[i];
	delete[] final_subchannel_collection;

	for (int i = 0; i < threads; ++i)
		delete[]  final_subchannel_pch_num[i];
	delete[] final_subchannel_pch_num;

	delete[] final_channel_pch_num;
	delete[] final_pch_num;
	delete[] final_subchannel_num;

	return;
}//END OF ASP