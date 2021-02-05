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

void channel_parallel_all(
	struct patch_object* patch,
	int available_threads,
	double xita_g,
	double Tcg,
	struct status_information* INF,
	bool verbose
) {


	int layer = 0;
	int init_flag = 0; // 0, temporal; not 0, fianl result


	init_channel_parallel(patch, INF);//make preparations before parallelization

	//start parallelization
	//the number of layers and iteration times is uncertain
	do
	{
		//local parameters
		double Tchm = Tcg;//set as the smallest unit
		int nlayer_potential_schemes = 0;//recond the number of potential schemes
		int maximum_nlayer_potential_schemes = 1;
		do {
			maximum_nlayer_potential_schemes = maximum_nlayer_potential_schemes + 2;
			Tchm = max(Tchm * (1 + xita_g), Tchm + Tcg);
		} while (Tchm <= INF->channel_num);
		// cout << maximum_nlayer_potential_schemes << endl;

		Tchm = Tcg;//reset as the smallest unit
		struct status_information* INFch = new struct status_information[maximum_nlayer_potential_schemes];//store INF for each potential schemes
		int terminate_layer_flag = 0;//to break out the iteration of layers

		//start parallelization for a layer
		do {

			//generate a potential scheme by channel unit partition and thread allocation
			channel_parallel(patch, Tchm, available_threads, INF);

			//cout informations
			if (verbose)
				cout << "Threads:" << available_threads << "\tlayer:" << layer << "\tTchm:" << fixed << setprecision(2) << Tchm << "\t nlayer_potential_schemes:" << nlayer_potential_schemes
				<< "\tESRch_layer:" << fixed << setprecision(2) << INF->ESRch << "\tTOV:" << fixed << setprecision(2) << INF->TOV << endl;

			//return to the original state of previous layers with 0
			init_flag = 0; // a temporal try with 0
			update_channel_parallel_within_layer(patch, layer, init_flag, INF);

			//store INF information
			INFch[nlayer_potential_schemes] = INF[0];

			//check for termination of this layer
			if (nlayer_potential_schemes > 1) {
				//see Fig.4 in the manuscript
				//it is essential to terminate the iteration processes as SR drops
				if (INFch[nlayer_potential_schemes].ESRch - INFch[nlayer_potential_schemes - 1].ESRch < 0)
					terminate_layer_flag = 1;
				else if (INF->rest_num == 0)
					terminate_layer_flag = 1;
			}

			//move to next potential scheme
			Tchm = max(Tchm * (1 + xita_g), Tchm + Tcg);
			nlayer_potential_schemes++;

		} while (terminate_layer_flag == 0);


		// select optimal result in terms of ESRch and TOV
		double max_ESRch = 0;
		double max_TOV = 0;
		int ioptimal_scheme = 0;
		//search max_ESRch and max_TOV
		for (int inx = 0; inx != nlayer_potential_schemes; inx++)
			if (INFch[inx].ESRch > max_ESRch)
				max_ESRch = INFch[inx].ESRch;
		for (int inx = 0; inx != nlayer_potential_schemes; inx++)
			if (fabs(INFch[inx].ESRch - max_ESRch) < 0.0001 && INFch[inx].TOV > max_TOV) {
				max_TOV = INFch[inx].TOV;
				ioptimal_scheme = inx;
			}
		delete[] INFch;

		// recompute Tchm and generate the optimal scheme
		Tchm = Tcg;
		for (int inx = 0; inx != ioptimal_scheme; inx++)
			Tchm = max(Tchm * (1 + xita_g), Tchm + Tcg);
		channel_parallel(patch, Tchm, available_threads, INF);
		init_flag = 1; //set as fianl results with 1
		update_channel_parallel_within_layer(patch, layer, init_flag, INF);

		//store informations for this thread 
		INF->TOV_layers[layer] = INF->TOV;
		INF->ESRch_layers[layer] = INF->ESRch; // set ESR to the layers

		//cout informations
		if (verbose)
			cout << "Threads:" << available_threads << "\tlayer:" << layer << "\tTchm:" << fixed << setprecision(2) << Tchm << "\t nlayer_potential_schemes:" << ioptimal_scheme
			<< "\tESRch_layer:" << fixed << setprecision(2) << INF->ESRch_layers[layer] << "\tTOV:" << fixed << setprecision(2) << INF->TOV << endl;

		//mover to next layer
		layer++;

	} while (INF->rest_num != 0);

	INF->LAYER = layer;

	//compute final ESRch
	INF->ESRch = 0;
	for (int inx = 0; inx != INF->LAYER; inx++)
		INF->ESRch += INF->TOV_layers[inx] / INF->ESRch_layers[inx];

	INF->ESRch = INF->channel_num / INF->ESRch;

	return;

}