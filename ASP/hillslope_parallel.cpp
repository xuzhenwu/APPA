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

void hillslope_parallel(struct patch_struct *patch, double Deviation_Index,int Decomposition_Index,  int threads, struct inf *INF,int out_flag) {

	//=======================================================================================================================
	//1. LOCAL VARS(JUST SKIP IF YOU DONT WANT TO READ ABOUT IT)
	//=======================================================================================================================

	const int infinite = 99999999;
	// while partitioning, num of subbasins will be more than the basin_num we want
	// thus we need more memory for storage and by using optimize method we finally got sufficient basin_num
	
	int itmp{};
	int break_flag{};

	// maximum runtime of sub-basins
	// the limit of accumulated runtime of a new basin
	double SM = INF->hillslope_num / (Decomposition_Index *threads) * (1 + Deviation_Index);

	//an array stores the sub-basins and patches in partition process
	const int sub_basin_times(20);//set as 20 maximum times of threads
	int **subbasin_grids_id = new int *[sub_basin_times * threads];
	for (int i = 0; i < sub_basin_times * threads; ++i)
	{
		subbasin_grids_id[i] = new int[INF->hillslope_num]{};
	}
	//coresponding sub-basin id and grids_num 
	int *subbasin_id = new int[sub_basin_times * threads]{};
	int *subbasin_grids_num = new int[sub_basin_times * threads]{};
	double *subbasin_hill_acc= new double[sub_basin_times * threads]{};

	//=======================================================================================================================
	//2. STEP 1:SUB-BASIN PARTITION
	//=======================================================================================================================
	// the ID of basin num 
	int basin_inx = 0;//start from 0 and controls the num of subbasins (before we optimize it)
	int rest_num = INF->patch_num;;//rest of num undistributed with their basin_inx

	do {
		//local parameters
		int satisfied_subbasin_hi_acc = 0;
		int satisfied_subbasin_outlet_ID = 0;
		int satisfied_deviation = infinite;

		//search for a satisfied sub-basin
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
		} 

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

		//to terminate the partition porcess
		if (out_flag != 0)
			rest_num = 0;

	} while (rest_num != 0);

	INF->subbasin_num = basin_inx;
	
	//reinit the sub-basins array
	//exclue channel grids in the list
	for (int inx = 0; inx < INF->subbasin_num; inx++) {
		subbasin_grids_num[inx] = 0;
	}
	for (int pch = 0; pch != INF->patch_num; pch++) {
	
		if (patch[pch].basin_inx > 0 && patch[pch].landID == 0) {
			subbasin_grids_id[patch[pch].basin_inx-1][(subbasin_grids_num[patch[pch].basin_inx])] = patch[pch].pch;//in current version we use pch to accrelate output
			subbasin_grids_num[patch[pch].basin_inx - 1]++;
		}
	}

	//=======================================================================================================================
	//3. STEP 2: THREAD ALLOCATION. allocate sub-basins into final_threads 
	//=======================================================================================================================
	
	//sub-basins in threads
	int **final_subbasin_id_in_threads = new int *[threads] {};
	for (int i = 0; i < threads; ++i)
	{
		final_subbasin_id_in_threads[i] = new int[threads*sub_basin_times]{};
	}
	//num of patches in a subbasin 
	int **final_subbasin_grids_num_in_threads = new int *[threads] {};
	for (int i = 0; i < threads; ++i)
	{
		final_subbasin_grids_num_in_threads[i] = new int[threads*sub_basin_times]{};
	}
	//num of patches in a subbasin 
	int **final_subbasin_hill_acc_in_threads = new int *[threads] {};
	for (int i = 0; i < threads; ++i)
	{
		final_subbasin_hill_acc_in_threads[i] = new int[threads*sub_basin_times]{};
	}

	//optimizations
	//num of patches in a basin
	int *final_grids_num_in_threads = new int [threads] {};
	//num of subbasins in a basin
	int *final_subbasin_num_in_threads = new int [threads] {};
	//num of subbasins in a basin
	double *final_hill_acc_in_threads = new double [threads] {};


	//final girds in thread
	int *final_grids_num_in_threads = new int [threads] {};

	//add basin_ID
	for (int inx = 0; inx < basin_inx; inx++)
		subbasin_id[inx] = inx;

	//local parameters
	int min_subsubbasin_grids_num{};
	int min_subbasin_pch_ID{};
	int max_flag{};
	int max_basin_ID{};
	int max_subbasin_ID_OUT{};
	int max_subbasin_grids_num{};
	int optimize_total{};
	int optimize_flag{};


	// put basins in the final basin assemble
	for (int inx = 0; inx < INF->subbasin_num; inx++) {

		int thread_inx= inx % threads;
		final_subbasin_id_in_threads[thread_inx][(final_subbasin_num_in_threads[thread_inx])] = subbasin_id[inx];//from 0
		final_subbasin_grids_num_in_threads[thread_inx][(final_subbasin_num_in_threads[thread_inx])] = subbasin_grids_num[inx];
		final_subbasin_hill_acc_in_threads[thread_inx][(final_subbasin_num_in_threads[thread_inx])] = subbasin_grids_num[inx];
		final_hill_acc_in_threads[thread_inx]+=subbasin_hill_acc[inx];
		final_subbasin_num_in_threads[thread_inx]++;
		final_grids_num_in_threads[thread_inx] += subbasin_grids_num[inx];
	}

	//main optimizing process to search for an ideal subbasin assemble to final basin
	do {

		//reinit as 0
		int maximum_hill_acc_in_threads_id = 0;
		double maximum_hill_acc_in_threads = 0;

		for (int inx = 0; inx != threads; inx++) {

			if (final_hill_acc_in_threads[inx] > maximum_hill_acc_in_threads) {
				maximum_hill_acc_in_threads_id = inx; 
				maximum_hill_acc_in_threads = final_hill_acc_in_threads[inx];
			}
		}
		for()





		//find it's smallest subbasin
		for (int iny = 0; iny != final_subbasin_num_in_threads[max_basin_ID]; iny++) {

			if (iny == 0) {

				min_subsubbasin_grids_num = final_subbasin_grids_num_in_threads[max_basin_ID][iny];
				max_subbasin_ID_OUT = iny;

			}

			if (final_subbasin_grids_num_in_threads[max_basin_ID][iny] < min_subsubbasin_grids_num) {

				min_subsubbasin_grids_num = final_subbasin_grids_num_in_threads[max_basin_ID][iny];
				max_subbasin_ID_OUT = iny;
			}

		}

		//check if we can optimize it 
		//by adding it to one of other basins
		int optimize_flag = 0;

		for (int inx = 0; inx != threads; inx++) {

			//only in other basins
			if (inx != max_basin_ID) {

				if (final_grids_num_in_threads[inx] + min_subsubbasin_grids_num < final_grids_num_in_threads[max_basin_ID]) {

					//remove smallest subbasins into this basin

					//an -1 as it's ID is from 0
					final_subbasin_id_in_threads[inx][(final_subbasin_num_in_threads[inx])] = final_subbasin_id_in_threads[max_basin_ID][max_subbasin_ID_OUT];
					final_subbasin_grids_num_in_threads[inx][(final_subbasin_num_in_threads[inx])] = final_subbasin_grids_num_in_threads[max_basin_ID][max_subbasin_ID_OUT];

					final_subbasin_id_in_threads[max_basin_ID][max_subbasin_ID_OUT] = 0;
					final_subbasin_grids_num_in_threads[max_basin_ID][max_subbasin_ID_OUT] = 0;

					//need to reorder the list of subbasins in this OUT basin
					for (int ainx = 0; ainx != final_subbasin_num_in_threads[max_basin_ID]; ainx++)
						for (int binx = ainx + 1; binx != final_subbasin_num_in_threads[max_basin_ID]; binx++) {

							if (final_subbasin_grids_num_in_threads[max_basin_ID][ainx] < final_subbasin_grids_num_in_threads[max_basin_ID][binx]) {

								//exchange value
								itmp = final_subbasin_id_in_threads[max_basin_ID][ainx];
								final_subbasin_id_in_threads[max_basin_ID][ainx] = final_subbasin_id_in_threads[max_basin_ID][binx];
								final_subbasin_id_in_threads[max_basin_ID][binx] = itmp;

								itmp = final_subbasin_grids_num_in_threads[max_basin_ID][ainx];
								final_subbasin_grids_num_in_threads[max_basin_ID][ainx] = final_subbasin_grids_num_in_threads[max_basin_ID][binx];
								final_subbasin_grids_num_in_threads[max_basin_ID][binx] = itmp;

							}
						}

					//then we can delete the value
					final_subbasin_num_in_threads[inx]++;
					final_subbasin_num_in_threads[max_basin_ID]--;
					final_grids_num_in_threads[inx] += min_subsubbasin_grids_num;
					final_grids_num_in_threads[max_basin_ID] -= min_subsubbasin_grids_num;

					optimize_flag = 1;
					optimize_total++;

					//cout << "round\t" << optimize_total << endl;

					//need to break 
					break;

				}
			}
		}


	} while (optimize_flag == 1);// till we got fine size of final basins or it can not be optimized


	//correlating final basins with final thread
	for (int pch = 0; pch != INF->patch_num; pch++) {

		//search
		break_flag = 0;
		for (int basin_inx = 0; basin_inx != threads; basin_inx++)
		{

			for (int sub_basin_inx = 0; sub_basin_inx != final_subbasin_num_in_threads[basin_inx]; sub_basin_inx++) {

				if (patch[pch].basin_inx == (final_subbasin_id_in_threads[basin_inx][sub_basin_inx] + 1)) {

					//start from 1
					//patch[pch].basin_inx = basin_inx + 1; do not update any more
					patch[pch].sthread = basin_inx + 1;
					
					final_grids_num_in_threads[basin_inx]++;

					//break in this for and next for
					break_flag = 1;
					break;

				}

			}
			//which means we has searched a good point
			if (break_flag == 1) break;

		}
		if (break_flag != 1) cout << "UNMSMCHED PSMCH WITH THEIR BASIN INX" << endl;
	}

	//CHECK DF
	double Final_Deviation = 0;
	for (basin_inx = 0; basin_inx != threads; basin_inx++) {
		if (final_grids_num_in_threads[basin_inx] > Final_Deviation)
			Final_Deviation = final_grids_num_in_threads[basin_inx];
	}
	Final_Deviation = Final_Deviation / (INF->hillslope_num / threads) - 1;


	//FREE MEMORY
	for (int i = 0; i < sub_basin_times * threads; ++i)
			delete [] subbasin_grids_id[i];
	delete[] subbasin_grids_id;

	delete[] subbasin_id;
	delete[] subbasin_grids_num;

	for (int i = 0; i < threads; ++i)
		delete[] final_subbasin_id_in_threads[i];
	delete[] final_subbasin_id_in_threads;
	
	for (int i = 0; i < threads; ++i)
		delete[]  final_subbasin_grids_num_in_threads[i];
	delete[] final_subbasin_grids_num_in_threads;

	delete[] final_grids_num_in_threads;
	delete[] final_grids_num_in_threads;
	delete[] final_subbasin_num_in_threads;

	INF->DF = Final_Deviation;
	INF->SR = 1/(1+Final_Deviation)*threads;

	return;
}//END OF ASP