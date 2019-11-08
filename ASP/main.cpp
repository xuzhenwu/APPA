//##########################################
//## MULTITHREAD PROCESSING FOR CELL GROUPS
//##########################################

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

//=======================================================================================================================
//1. DEFINE TUNABLE PARAMETERS
//=======================================================================================================================

/*
//file location
char flowtable_in[121]("I://Data//CHESS//IN//xf_ws//geo//xf_ws_flow_table_D8.dat");// flow table in
char patch_in[121]("I://Data//CHESS//IN//xf_ws//geo//xf_ws.patch");// .patch file in, patchID
char sthread_out[121]("I://Data//CHESS//IN//xf_ws//geo//xf_ws.sthread");//sub-basin thread
char cthread_out[121]("I://Data//CHESS//IN//xf_ws//geo//xf_ws.cthread");//channel thread
*/

//char geo_files[121]("I://xu//CHESS_PARALLEL//Data//xf_ws//geo//xf_ws");
//char geo_files[121]("I://xu//CHESS_PARALLEL//Data//lc//geo//lc");
char geo_files[121]("Z://PARA//xf_ws//geo//xf_ws");
double  SDI(8);//Sub-basin division index, it will increase more sub-basins when it excels 1

double  Resolution(200);//m

//thread parameters
const int thread_num = 17;
int threads[thread_num] = {1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32 };

double  MAX_DT = 0.4;//maximum DT
double  RDT = 0.05;// increment for next DT

double  PDT(0.2);//used in channel partition
double  SPCU(0.25);//smallest_portion_channel_unit
double  SR_GOAL_PER(1);


//other parameters (less important)
double  Stream_Boun_Area(0.5);//km^2, this is for refined landID of a patch when needed
const int max_out_inf(100000);
struct inf *OUT_INF = new struct inf [max_out_inf];
struct inf *OUT = new struct inf;

int layer;


using namespace std;

int main()
{
	int patch_num;
	int *ppatch_num=new int;
	char flowtable_in[121]{}, patch_in[121]{}, sthread_out[121]{}, cthread_out[121]{};
	strcpy_s(flowtable_in, geo_files);
	strcat_s(flowtable_in, "_flow_table_D8.dat");

	strcpy_s(sthread_out, geo_files);
	strcat_s(sthread_out, ".sthread");
	strcpy_s(cthread_out, geo_files);
	strcat_s(cthread_out, ".cthread");
	cout << SDI << endl;
	// allow for mutiple thread chioces
	for (int thread_inx = 0; thread_inx !=17; thread_inx++) {
		//for (int thread_inx = 0; thread_inx != 8; thread_inx++) {
			//for (int thread_inx = 8; thread_inx != 17; thread_inx++) {

		//print into screen
		//cout << "\tTHREAD " << threads[thread_inx] <<"\t FLOW TABLE"<< endl;

		//Build and Read flow table 
		struct patch_struct *patch = NULL;
		patch= read_flow_table(patch, flowtable_in, OUT);
		patch_num = OUT->patch_num;

		///*
		//=======================================================================================================================
		// SUB-BASIN PARTION AND THREAD COMPOSITION
		//=======================================================================================================================
		int sub_try_times = 0; 
		double DT = 0;//the first value of deviation_index
		do
		{

			init_hillslope_parallel(patch, OUT);

			//sub-bain partition
			hillslope_parallel(patch, DT, SDI, threads[thread_inx], OUT,0);

			cout << "\t\t\t " << DT << "\t"  << fixed << setprecision(2) << OUT->SR << "\t" << fixed << setprecision(2) << OUT->DF
				<< "\t" << OUT->subbasin_num  << endl;

			//init
			


			DT += RDT;
			OUT_INF[sub_try_times] = OUT[0];
			sub_try_times++;

		} while (DT < MAX_DT);
		
		
		//select optimal result by goals
		double max_sr = 0;
		int best_OUT_inx = 0;
		//[SR first]in default
		for (int inx = 0; inx != sub_try_times; inx++) {
			if (OUT_INF[inx].SR > max_sr) {
				max_sr = OUT_INF[inx].SR;
				best_OUT_inx = inx;
			}
		}

		DT = best_OUT_inx * RDT;

		//sub-bain partition
		hillslope_parallel(patch, DT, SDI, threads[thread_inx], OUT,0);

		cout << "\t\t\t " << DT << "\t" << fixed << setprecision(2) << OUT->SR << "\t" << fixed << setprecision(2) << OUT->DF
			<< "\t" << OUT->subbasin_num << endl;

		OUT->HSR = OUT->SR;

		//WRITE STHREAD FILE 
		print_threads_map(patch, geo_files, sthread_out, patch_num, threads[thread_inx], 1);





		//=======================================================================================================================
		// CHANNEL PARTION AND THREAD COMPOSITION
		//=======================================================================================================================
		double SR_GOAL = threads[thread_inx] * SR_GOAL_PER;		
		layer = 0;

		do
		{
			//init
			DT = 0;
			int layer_try_times = 0;

			//iteration in each layer
			int out_flag = 0;
			do{

				//computing all informations
				channel_parallel(patch, DT, SDI, threads[thread_inx], patch_num, Resolution, SPCU, OUT,0);
				
				//cout informations
				//cout << "\t " << layer <<"\t"<<layer_try_times<< "\t" << fixed << setprecision(2) << DT << "\t" << fixed << setprecision(2) << OUT->SR
					//<< "\t" <<OUT->CUN << "\t" << OUT->PN << "\t" << OUT->RE << "\t" << OUT->CN << endl;
				
				//init as in original state
				init_channel_parallel(patch, patch_num, layer, 0);

				//store out information
				OUT_INF[layer_try_times]=OUT[0];

				//move to next layer try
				DT = (1+DT)*(1+RDT)-1;
				
				if (layer_try_times > 1) {
					if (OUT_INF[layer_try_times].SR - OUT_INF[layer_try_times - 1].SR < 0)
						out_flag = 1;
					else if (OUT->RE == 0)
						out_flag = 1;
				}

				layer_try_times++;
			} while (out_flag==0);//


			//select optimal result by goals
			double max_sr = 0;
			double max_pn = 0;
			int best_OUT_inx=0;
			//compute max sr
			for (int inx = 0; inx != layer_try_times; inx++) {
				
				if (OUT_INF[inx].SR > max_sr){
					max_sr = OUT_INF[inx].SR;
					best_OUT_inx = inx;
				}
			}
			//[PN first]
			if (max_sr >= SR_GOAL) {
				for (int inx = 0; inx != layer_try_times; inx++) {
					//search_max_sr with max_pn
					if (OUT_INF[inx].SR >= SR_GOAL && OUT_INF[inx].PN > max_pn) {//excels to SR_GOAL with max PN
						max_pn = OUT_INF[inx].PN;
						best_OUT_inx = inx;
					}
				}
			}
			//[SR first]
			else {
			
				for (int inx = 0; inx != layer_try_times; inx++) {
					if (fabs(OUT_INF[inx].SR-max_sr)<0.0001 && OUT_INF[inx].PN > max_pn) {//equals to MAX_SR with max PN
						max_pn = OUT_INF[inx].PN;
						best_OUT_inx = inx;
					}
				}
			
			}

			//renew layer result
			DT = 0;
			for (int inx = 0; inx != best_OUT_inx;inx++)
				DT = (1 + DT)*(1 + RDT)-1;

			channel_parallel(patch, DT, SDI, threads[thread_inx], patch_num, Resolution, SPCU, OUT,0);

			//cout informations
			cout << "\t " << layer << "\t" << best_OUT_inx <<"\t"<< fixed << setprecision(2) << DT << "\t" << fixed << setprecision(2) << OUT->SR
				<< "\t" << OUT->CUN << "\t" << OUT->PN << "\t" << OUT->RE << "\t" << OUT->CN << endl;

			//as static
			init_channel_parallel(patch, patch_num, layer, 1);
			OUT->PN_CH[layer] = OUT->PN;
			OUT->SR_CH[layer] = OUT->SR;
			layer++;

			if (layer == 1) {
				int a = 0;
				//break;
			}


		} while (OUT->RE!=0);
		OUT->LAYER = layer;

		//print channel thread map
		//print_threads_map(patch, geo_files, cthread_out, patch_num, threads[thread_inx], 2);

		//compute final CSR
		for (int inx = 0; inx != OUT->LAYER; inx++) {
		
			OUT->CSR += OUT->PN_CH[inx] / OUT->SR_CH[inx];
		}
		OUT->CSR = OUT->CN / OUT->CSR;

		//print results
		
		OUT->SR = patch_num / ((patch_num - OUT->CN) / OUT->HSR + OUT->CN / OUT->CSR);
		
		cout << "\tTHREAD " << threads[thread_inx] << "\t End with " <<  fixed << setprecision(2) << OUT->SR<<  "\t"<<
			fixed << setprecision(2) << OUT->HSR << "\t" << patch_num- OUT->CN << "\t" <<
			fixed << setprecision(2) << OUT->CSR << "\t" << OUT->CN << "\t\n" << endl;

		delete patch;

	}

	

	cout << "\nOUTPUT IS AVAILABLE NOW...." << endl;
	getchar();

	return 0;
}