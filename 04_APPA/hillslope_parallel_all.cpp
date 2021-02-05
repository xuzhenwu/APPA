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

void hillslope_parallel_all(
	struct patch_object *patch, 
	int available_threads,
	double xita_s_range[3], 
	int D_range[3],
	struct status_information *INF,
	bool verbose 
) 
{

	//=======================================================================================================================
	// 3-1 Hillslope parallelization
	//=======================================================================================================================

	int nD = (D_range[1] - D_range[0]) / D_range[2] + 1;//number of D
	int n_xita_s = floor((xita_s_range[1] - xita_s_range[0]) / xita_s_range[2]) + 1;//number of xita_s
	double max_ESRhs= 0;
	double optimal_D, optimal_xita_s;

	//search for potential schemes from D and xita_s ranges
	for (int iD = 0; iD < nD; iD++)
		for (int i_xita_s = 0; i_xita_s < n_xita_s; i_xita_s++) {

			//local parameters
			double xita_s = xita_s_range[0] + i_xita_s * xita_s_range[2];
			double D = D_range[0] + iD * D_range[2];

			//make preparations before parallelization
			init_hillslope_parallel(patch, INF);
			//generate a potential scheme by sub-bain partition and thread allocation
			hillslope_parallel(patch, xita_s, D, available_threads, INF);

			//print informations to screens
			if (verbose != 0)
				cout << "Threads:" << available_threads << "\txita_s:" << fixed << setprecision(2) << xita_s << "\tD:" << D
				<< "\tESRhs:" << fixed << setprecision(2) << INF->ESRhs<< "\tm:" << fixed << setprecision(2) << INF->subbasin_num << endl;

			//select optimal result in terms of ESRhs
			if (INF->ESRhs > max_ESRhs) {
				max_ESRhs = INF->ESRhs;
				optimal_D = D;
				optimal_xita_s = xita_s;
			}
		}

	// generate the optimal scheme
	init_hillslope_parallel(patch, INF);
	hillslope_parallel(patch, optimal_xita_s, optimal_D, available_threads, INF);

	//print informations to screens
	if (verbose)
		cout << "Threads:" << available_threads << "\txita_s:" << fixed << setprecision(2) << optimal_xita_s << "\tD:" << optimal_D
		<< "\tESRhs:" << fixed << setprecision(2) << INF->ESRhs<< "\tm:" << fixed << setprecision(2) << INF->subbasin_num << endl;

	return;

}