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


void hillslope_parallel(struct patch_object *patch, double xita_s, int D,  int threads, struct status_information *INF){

	//=======================================================================================================================
	//1. ENVIRONMENT SETTING
	//=======================================================================================================================
	// Eq. 5 and 6 
	double Thsn = INF->hillslope_num / threads;
	double SM = Thsn / D  * (1 + xita_s); // maximum computation time demanded for partitioned sub-basins

	//an array stores the sub-basins and patches in partition process
	int sub_basin_times = D*2;//set as 2 maximum times of threads*D, it may exceed the limit while using a too-large value of D for a small basin
	int **subbasin_grids_id = new int *[sub_basin_times * threads];

	int max_subbasin_size = SM * 4;

	for (int i = 0; i < sub_basin_times * threads; i++)
	{
		subbasin_grids_id[i] = new int[max_subbasin_size]{};//use double memory for avoiding break ups 
	}
	//coresponding sub-basin id and grids_num
	int *subbasin_id = new int[sub_basin_times * threads]{};
	int *subbasin_grids_num = new int[sub_basin_times * threads]{};
	double *subbasin_hill_acc = new double[sub_basin_times * threads]{}; 
	
	const int infinite = 99999999;

	//=======================================================================================================================
	//2. STEP 1:SUB-BASIN PARTITION
	//=======================================================================================================================
	// the ID of basin num 
	int basin_inx = 0;//start from 0 and controls the num of subbasins (before we optimize it)
	int rest_num = INF->patch_num;;//number of the rest grids unallocated with their basin_inx

	do {
		//local parameters
		double satisfied_subbasin_hi_acc = 0;
		int satisfied_subbasin_outlet_ID = 0;
		double satisfied_deviation = infinite;

		//search for a satisfied sub-basinstart sorting in terms of elevation
		for (int pch = 0; pch != INF->patch_num; pch++) {

			//stream grid that are remained 
			if (patch[pch].landID ==1 && patch[pch].basin_inx == 0)
			{
				//find best ID
				if (fabs(patch[pch].re_hill_acc - SM)< satisfied_deviation && patch[pch].re_hill_acc - SM<0) {

						satisfied_subbasin_hi_acc = patch[pch].re_hill_acc;
						satisfied_subbasin_outlet_ID = pch;
						satisfied_deviation = fabs(patch[pch].re_hill_acc - SM);
				}
			}
		} // an outlet with tgh lower than SM
		if (satisfied_subbasin_outlet_ID == 0) {
			for (int pch = 0; pch != INF->patch_num; pch++)
			if (patch[pch].landID == 1 && patch[pch].basin_inx == 0)
			{
				//find best ID
				if (fabs(patch[pch].re_hill_acc - SM)< satisfied_deviation && patch[pch].re_hill_acc - SM>0) {

					satisfied_subbasin_hi_acc = patch[pch].re_hill_acc;
					satisfied_subbasin_outlet_ID = pch;
					satisfied_deviation = fabs(patch[pch].re_hill_acc - SM);
				
				}
			}
			if(satisfied_subbasin_hi_acc >= max_subbasin_size){
				cout << "an outlet with Tg,s higher than SM\t" << SM << "\t Partition\t" << satisfied_subbasin_hi_acc << endl;
				cout << "this is limited by the extendence of channel routing" << endl;
				cout << "please switch to a small value of 'D' and restart the program !";
				exit(0);
			}
		}// an outlet with tgh higher than SM

		//distribute basin_inx and update re_ch_acc_area
		patch[satisfied_subbasin_outlet_ID].basin_inx = basin_inx + 1;//start from 1
		subbasin_grids_id[basin_inx][(subbasin_grids_num[basin_inx])] = patch[satisfied_subbasin_outlet_ID].patchID;
		subbasin_grids_num[basin_inx]++;
		subbasin_hill_acc[basin_inx] = satisfied_subbasin_hi_acc;

		//searching for upslope undefined patches
		for (int up_pch_inx = satisfied_subbasin_outlet_ID - 1; up_pch_inx >= 0; up_pch_inx--) {

			// only when it's unlabelled
			if (patch[up_pch_inx].basin_inx == 0) {
				//check it's neigh
				for (int neigh = 0; neigh != patch[up_pch_inx].neigh_num; neigh++) {

					//check the list now in this basins
					for (int basin_pch_inx = 0; basin_pch_inx != subbasin_grids_num[basin_inx]; basin_pch_inx++) {

						//when it's neigh is in the basin list
						if (patch[up_pch_inx].neigh_ID[neigh] == subbasin_grids_id[basin_inx][basin_pch_inx]) {

							//add label
							patch[up_pch_inx].basin_inx = basin_inx + 1;//start from 1
							patch[up_pch_inx].re_hill_acc = 0;

							//add it's ID in the basin list when it is hillslope grids
							subbasin_grids_id[basin_inx][(subbasin_grids_num[basin_inx])] = patch[up_pch_inx].patchID;
							subbasin_grids_num[basin_inx]++;

							break;
						}
					}

				}
			}

		}//end of seaching it's upslope area

		//after divided
		rest_num -= subbasin_grids_num[basin_inx];
		basin_inx++;//move for next sub-basin

		//renew downstream re_hill_acc
		if (rest_num != 0) {
		
			int start_pch = satisfied_subbasin_outlet_ID;
			do {
				patch[start_pch].re_hill_acc -= satisfied_subbasin_hi_acc;
				if(patch[start_pch].neigh_num>0)
				   start_pch = patch[start_pch].neigh_pch[0];
				else break;//when there is no neighbor grids
			} while (start_pch <= INF->patch_num);
		}


		if (basin_inx > sub_basin_times*threads)
			printf("ERROR: the number of basin_inx is larger than the limit, %d > %d * %d nthreads", basin_inx, sub_basin_times, threads);

	} while (rest_num != 0);

	INF->subbasin_num = basin_inx;
	
	//store subbasin_id
	for (int inx = 0; inx < basin_inx; inx++)
		subbasin_id[inx] = inx;

	//reinit the sub-basins array
	//exclue channel grids in the list
	for (int inx = 0; inx < INF->subbasin_num; inx++) {
		subbasin_grids_num[inx] = 0;
	}
	for (int pch = 0; pch != INF->patch_num; pch++) {
	
		if (patch[pch].basin_inx > 0 && patch[pch].landID == 0) {
			subbasin_grids_id[patch[pch].basin_inx-1][(subbasin_grids_num[patch[pch].basin_inx-1])] = patch[pch].pch;//in current version we use pch to accrelate output
			subbasin_grids_num[patch[pch].basin_inx - 1]++;
		}
	}

	//=======================================================================================================================
	//3. STEP 2: THREAD ALLOCATION WITH LOAD BALANCE  
	//=======================================================================================================================

	//sub-basins id in threads
	int **final_subbasin_id_in_threads = new int *[threads] {};
	for (int i = 0; i < threads; ++i)
	{
		final_subbasin_id_in_threads[i] = new int[threads*sub_basin_times]{};
	}
	//sub-basins hill_acc in threads 
	double **final_subbasin_hill_acc_in_threads = new double *[threads] {};
	for (int i = 0; i < threads; ++i)
	{
		final_subbasin_hill_acc_in_threads[i] = new double[threads*sub_basin_times]{};
	}

	//sum up for threads
	//num of subbasins in a basin
	int *final_subbasin_num_in_threads = new int [threads] {};
	//num of subbasins in a basin
	double *final_hill_acc_in_threads = new double [threads] {};

	int optimize_flag{};


	// put basins in the final basin assemble
	for (int inx = 0; inx < INF->subbasin_num; inx++) {

		int thread_inx= inx % threads;
		final_subbasin_id_in_threads[thread_inx][(final_subbasin_num_in_threads[thread_inx])] = subbasin_id[inx];//from 0
		final_subbasin_hill_acc_in_threads[thread_inx][(final_subbasin_num_in_threads[thread_inx])] = subbasin_grids_num[inx];
		final_hill_acc_in_threads[thread_inx]+=subbasin_hill_acc[inx];
		final_subbasin_num_in_threads[thread_inx]++;
	}

	//main optimizing process to search for an ideal subbasin assemble to final threads
	do {

		//find maximum thread
		int maximum_hill_acc_in_threads_id = 0;
		double maximum_hill_acc_in_threads = 0;
		int minimum_hill_acc_in_threads_id = 0;
		double minimum_hill_acc_in_threads = infinite;
		for (int inx = 0; inx != threads; inx++) {

			if (final_hill_acc_in_threads[inx] > maximum_hill_acc_in_threads) {
				maximum_hill_acc_in_threads_id = inx; 
				maximum_hill_acc_in_threads = final_hill_acc_in_threads[inx];
			}
			if (final_hill_acc_in_threads[inx] < minimum_hill_acc_in_threads) {
				minimum_hill_acc_in_threads_id = inx;
				minimum_hill_acc_in_threads = final_hill_acc_in_threads[inx];
			}
		}

		//find minimum sub-basin
		int minimum_subbasin_id = 0;
		double minimum_subbasin_hill_acc = infinite;
		for (int inx = 0; inx != final_subbasin_num_in_threads[maximum_hill_acc_in_threads_id]; inx++) {
		
			if (final_subbasin_hill_acc_in_threads[maximum_hill_acc_in_threads_id][inx] < minimum_subbasin_hill_acc) {
				minimum_subbasin_hill_acc = final_subbasin_hill_acc_in_threads[maximum_hill_acc_in_threads_id][inx];
				minimum_subbasin_id= inx;	
			}
		}

		//try to resdistribute to the minimum threads
		if (maximum_hill_acc_in_threads > minimum_hill_acc_in_threads + minimum_subbasin_hill_acc) {

			optimize_flag = 1;

			//move maximum thread, minimum subbasin to minimum thread
			//for minimum thread
			//add a new sub-basin
			final_subbasin_id_in_threads[minimum_hill_acc_in_threads_id][(final_subbasin_num_in_threads[minimum_hill_acc_in_threads_id])] =
				final_subbasin_id_in_threads[maximum_hill_acc_in_threads_id][minimum_subbasin_id];//move id of subbasin
			final_subbasin_hill_acc_in_threads[minimum_hill_acc_in_threads_id][(final_subbasin_num_in_threads[minimum_hill_acc_in_threads_id])]
				= minimum_subbasin_hill_acc;
			final_hill_acc_in_threads[minimum_hill_acc_in_threads_id] = minimum_hill_acc_in_threads + minimum_subbasin_hill_acc;//update hill_acc
			final_subbasin_num_in_threads[minimum_hill_acc_in_threads_id]++;//++ subbasin num 

			//for maximum thread
			//move id and hill_acc forward
			for (int inx = minimum_subbasin_id; inx != final_subbasin_num_in_threads[maximum_hill_acc_in_threads_id]; inx++) {
				final_subbasin_id_in_threads[maximum_hill_acc_in_threads_id][inx] = final_subbasin_id_in_threads[maximum_hill_acc_in_threads_id][inx + 1];
				final_subbasin_hill_acc_in_threads[maximum_hill_acc_in_threads_id][inx] = final_subbasin_hill_acc_in_threads[maximum_hill_acc_in_threads_id][inx + 1];
			}
			final_hill_acc_in_threads[maximum_hill_acc_in_threads_id] = maximum_hill_acc_in_threads - minimum_subbasin_hill_acc;
			final_subbasin_num_in_threads[maximum_hill_acc_in_threads_id]--;//++ subbasin num 

		}
		else
			optimize_flag = 0;


	} while (optimize_flag == 1);// till we got fine size of final basins or it can not be optimized

	//correlating final basins with final thread
	int partition_num = 0;
	for (int inx = 0; inx != threads; inx++) {

		for (int iny = 0; iny != final_subbasin_num_in_threads[inx]; iny++) {
		
			for (int inz = 0; inz != subbasin_grids_num[(final_subbasin_id_in_threads[inx][iny])]; inz++) {
			  
				int subbasin_id = final_subbasin_id_in_threads[inx][iny];
				patch[subbasin_grids_id[subbasin_id][inz]].sthread = inx+1;//start from 1
				partition_num++;
			}
		
		}
	}
	if (partition_num != INF->hillslope_num)
		cout << "partition_num do not equals to hillslope_num" << endl;


	//compute hsr
	double maximum_hill_acc = 0;
	double maximum_hill_acc_in_threads=0;
	for (int pch = 0; pch != INF->patch_num; pch++) {
		if (patch[pch].hill_acc > maximum_hill_acc) {
			maximum_hill_acc = patch[pch].hill_acc;
		}
	}
	for (int inx = 0; inx != threads; inx++) {
		if (final_hill_acc_in_threads[inx] >  maximum_hill_acc_in_threads)
			maximum_hill_acc_in_threads = final_hill_acc_in_threads[inx];
	}

	//Free MEMORY
	//local parameters in partition
	for (int i = 0; i < sub_basin_times * threads;i++)
			delete [] subbasin_grids_id[i];
	delete[] subbasin_grids_id;
	delete[] subbasin_id;
	delete[] subbasin_grids_num;
	delete[] subbasin_hill_acc;

	//local parameters in thread allocation
	for (int i = 0; i < threads; i++)
		delete[] final_subbasin_id_in_threads[i];
	delete[] final_subbasin_id_in_threads;
	for (int i = 0; i < threads; i++)
		delete [] final_subbasin_hill_acc_in_threads[i];
	delete[] final_subbasin_hill_acc_in_threads;

	delete [] final_subbasin_num_in_threads;
	delete [] final_hill_acc_in_threads;

	// status infomation
	INF->ESRhs = maximum_hill_acc / maximum_hill_acc_in_threads;
	INF->TMSRhs = threads;

	return;

}//end of hillslope parallelization