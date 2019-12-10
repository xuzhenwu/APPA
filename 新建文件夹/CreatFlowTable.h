#ifndef   CreatFlowTable_H_
#define   CreatFlowTable_H_


#include <iostream>
#include <iomanip>
#include <fstream>
#include "blender.h"
#include "CreatFlowTable.h"
#include "math.h"

using namespace std;

void CreatFlowTable::initialize() {
	/*double    *dem = { nullptr }, *slope = { nullptr }, *flna = { nullptr };
	int       *patch = { nullptr }, *stream = { nullptr }, *roads = { nullptr };*/
	flow_struct	*flow_table = { nullptr };
}

void CreatFlowTable::add_roads(flow_struct *flow_table, int num_patches, double cell)
{
	// local variable declarations
	int inx, neigh, pch, fnd;
	adj_struct *aptr;
	//ID_struct;

	// calculate gamma for each neighbour
	for (pch = 1; pch <= num_patches; pch++) {
		aptr = flow_table[pch].adj_list;
		if (flow_table[pch].land == 2) {
			fnd = find_stream(flow_table, pch, &inx);

			if (fnd == 0) exit(1);

			if (flow_table[inx].land == 1) {
				flow_table[pch].stream_ID.patch = flow_table[inx].patchID
					;
				flow_table[pch].stream_inx = inx;
				neigh = flow_table[pch].num_adjacent + 1;
			}
		}
	}

	return;
}

void CreatFlowTable::adjust_pit(flow_struct *flow_table, int curr, int edge_inx, double edge_elev, double cell)
{
	long int  j;
	double  total_perimeter;
	double  xrun, yrun;
	double  rise;
	adj_struct *aptr;

	aptr = flow_table[curr].adj_list;
	total_perimeter = aptr->perimeter;

	for (j = 1; j < flow_table[curr].num_adjacent; j++) {
		aptr = aptr->next;
		total_perimeter += aptr->perimeter;
	}

	aptr->next = new adj_struct;

	if (!aptr->next) {
		std::cout << "\n Not enough memory" << std::endl;
		exit(0);
	}

	aptr = aptr->next;
	aptr->patchID = flow_table[edge_inx].patchID;

	aptr->perimeter = total_perimeter;
	aptr->z = (double)edge_elev;
	aptr->inx = edge_inx;

	xrun = (flow_table[curr].x - flow_table[edge_inx].x);
	xrun = xrun * xrun;
	yrun = (flow_table[curr].y - flow_table[edge_inx].y);
	yrun = yrun * yrun;
	rise = flow_table[curr].z - flow_table[edge_inx].z;

	if ((yrun + xrun) <= 0.000000001)
		aptr->slope = rise;
	else
		//aptr->slope = ( double )(rise / (sqrt(xrun+yrun) * (cell) ) );
		aptr->slope = (double)(rise / sqrt(xrun + yrun));

	aptr->gamma = 1.0;
	flow_table[curr].total_gamma = (double)(aptr->perimeter *aptr->slope / (cell));
	flow_table[curr].z = (double)(edge_elev + MIN_RISE);
	flow_table[curr].num_adjacent += 1;

	return;
}

void build_flow_table(flow_struct *flow_table, double *dem, double *slope, int *patch, int *stream, int *roads, int *subbasin, double *flna, int maxr, int maxc, int f_flag, int sc_flag, double cell, double xllcorner, double yllcorner,
	int *patch_pch, int *rcboun, int thread_inx)
{
	// local variable declarations
	int inx{}, iny{};
	int r{}, c{}, pch{};
	double tmp_int{};
	double tmp_dou{};


	printf("\nSTART BUILDING FLOW TABLE\n");

	for (r = 0; r < maxr; r++) {
		for (c = 0; c < maxc; c++) {
			inx = r * maxc + c;

			if (inx >= rcboun[thread_inx] && inx < rcboun[thread_inx + 1]) {
				// TO: ignore areas outside the basin


				if (patch[inx] > 0) {

					//PCH STANDS FOR THE REAL ORDER NUM IN flow_table lists
					pch = patch_pch[inx];
					


					//pch start from 1
					flow_table[pch].patchID = patch[inx];

					flow_table[pch].subbasin_id = subbasin[inx];
					//TO: test for comparing algorithms acc_area here from one point
					/*if (flow_table[pch].patchID == 3400
						|| flow_table[pch].patchID == 11399
						|| flow_table[pch].patchID == 1410
						|| flow_table[pch].patchID == 10943
						){
						flow_table[pch].area = 1;
					}
					*/
					//INITIAL AREA AS 1
					flow_table[pch].area = 1;

					//INITIAL X,Y,Z
					flow_table[pch].x += (double)((1.0 * c)*cell + xllcorner);
					flow_table[pch].y += (double)(cell*(maxr - r - 1) + yllcorner);
					flow_table[pch].z += (double)(dem[inx]);


					if (sc_flag == 1)
						flow_table[pch].internal_slope += (double)(1.0 * slope[inx] * DtoR);

					//INITIAL LAND FLAG
					if (roads[inx] >= 1)
						flow_table[pch].land = 2;
					if (stream[inx] >= 1)
						flow_table[pch].land = 1;


					if (f_flag)
						flow_table[pch].flna += (double)(1.0 * flna[inx]);
					else
						flow_table[pch].flna = 0.0;

					//BUILD NEIGHBOUR LISTS
					flow_table[pch].num_adjacent += check_neighbours(r, c, patch,patch_pch, dem, stream,subbasin,
						&flow_table[pch],
						flow_table[pch].num_adjacent,
						maxr, maxc, sc_flag, cell);
				}//END OF A PATCH

			}
		}
	}
	cout << "END OF BUILDING FLOW TABLE of thread  " << thread_inx << " PATCHES" << endl;

	return;
}

int check_neighbours(int er, int ec, int *patch, int *patch_pch, double *dem, int *stream,int *subbasin, flow_struct *flow_entry, int num_adj, int maxr, int maxc, int sc_flag, double cell)
{
	/* local function declarations */
	//adj_struct  *check_list(int, int, adj_struct *);
	//int check_list_min();
	int r, c;
	int p_neigh;
	int stream_neigh;
	int pch_neigh;
	int new_adj = 0;
	int subbasin_neigh = 0;
	double ele_neigh = 0.;

	// for now streams don't point anywhere 			
	for (r = -1; r <= 1; r++) {
		for (c = -1; c <= 1; c++) {
			// don't look at neighbours beyond the edge 
			if ((er + r < maxr) && (er + r >= 0) && (ec + c < maxc) && (ec + c >= 0)) {

				// is the neighbour a different patch  or is it outside the basin - in which case we ignore it */
				// also, for stream pixels, ignore non-stream neighbours */
				p_neigh = patch[(er + r)*maxc + (ec + c)];
				stream_neigh = stream[(er + r)*maxc + (ec + c)];
				ele_neigh = dem[(er + r)*maxc + (ec + c)];
				pch_neigh = patch_pch[(er + r)*maxc + (ec + c)];
				subbasin_neigh = subbasin[(er + r)*maxc + (ec + c)];


				//Start is neigh
				if ((p_neigh != flow_entry->patchID) && (p_neigh > 0)) {

					/* create a list of neighbours if it does not exist already */
					if ((num_adj == 0)) {
						new_adj = 1;
						num_adj = 1;
						flow_entry->adj_list = new adj_struct;

						if (!flow_entry->adj_list)
						{
							std::cout << "\nMemory Allocation Failed for " << flow_entry->patchID << std::endl;
							exit(1);
						}

						flow_entry->adj_list->patchID = p_neigh;
						flow_entry->adj_list->z = ele_neigh;
						flow_entry->adj_list->inx = pch_neigh;
						flow_entry->adj_list->subbasin_id = subbasin_neigh;
						flow_entry->adj_list->perimeter = 0.0;
						flow_entry->adj_list->next = NULL;
						flow_entry->adj_ptr = flow_entry->adj_list;
					}
					else {
						// search list for other entries for this neighbour 
						flow_entry->adj_ptr = flow_entry->adj_list;
						flow_entry->adj_ptr = (struct adj_struct *)check_list(
							p_neigh, num_adj, flow_entry->adj_ptr);

						// is this the first instance of the neighbour
						if ((flow_entry->adj_ptr->patchID != p_neigh) && (p_neigh > 0))

						{
							// land processing 
							if (flow_entry->land != 1) { //&& (ele_neigh < flow_entry->z )
								flow_entry->adj_ptr->next = new adj_struct;
								if (!flow_entry->adj_ptr->next) {
									std::cout << "\nMemory Allocation Failed for " << flow_entry->patchID << std::endl;
									exit(1);
								}

								//printf_s("\nMemory Allocation Failed for %d %d", p_neigh,flow_entry->patchID);
								//printf_s("\nMemory Allocation Failed for %f %f", ele_neigh,flow_entry->z);
								//getchar();
								flow_entry->adj_ptr = flow_entry->adj_ptr->next;
								new_adj += 1;
								num_adj += 1;
								flow_entry->adj_ptr->patchID = p_neigh;
								flow_entry->adj_ptr->z = ele_neigh;
								flow_entry->adj_ptr->inx = pch_neigh;
								flow_entry->adj_list->subbasin_id = subbasin_neigh;
								flow_entry->adj_ptr->perimeter = 0.0;
								flow_entry->adj_ptr->next = NULL;
							} // end land processing 
							else { 		// stream processing
								if (stream_neigh != -9999) {  //exclude non-river type cells

									flow_entry->adj_ptr->next = new adj_struct;
									if (!flow_entry->adj_ptr->next) {
										std::cout << "\nMemory Allocation Failed for " << flow_entry->patchID << std::endl;
										exit(1);
									}

									flow_entry->adj_ptr = flow_entry->adj_ptr->next;

									new_adj += 1;
									num_adj += 1;
									flow_entry->adj_ptr->patchID = p_neigh;
									flow_entry->adj_ptr->z = ele_neigh;
									flow_entry->adj_ptr->inx = pch_neigh;
									flow_entry->adj_list->subbasin_id = subbasin_neigh;
									flow_entry->adj_ptr->perimeter = 0.0;
									flow_entry->adj_ptr->next = NULL;
								}
							}
						} // end new entry
					} //end of creating a list of neighbours

					// add perimeter length
					if (abs(c) + abs(r) == 1)
						flow_entry->adj_ptr->perimeter = 0.5;
					else
						flow_entry->adj_ptr->perimeter = 0.5 / sqrt(2.0);
				} // end is neigh
			} // end edges if
		} // end col 
	} // end row 

	flow_entry->adj_ptr = flow_entry->adj_list;

	return new_adj;
}

void CreatFlowTable::compute_drainage_density(flow_struct *flow_table, int num_patches, double cell)
{
	// local fuction declarations
	//struct ID_struct  sort_flow_table();
	//int	find_patch();

	// local variable declarations
	int inx;
	int neigh;
	int pch;
	int n_adjacent;
	double run, xrun, yrun;
	double drainage_density, total_area;
	struct adj_struct *aptr;
	struct ID_struct;

	total_area = 0.0;
	// send area to  each upslope neighbour for each patch in the sort list

	drainage_density = 0.0;
	n_adjacent = 0;

	for (pch = 1; pch <= num_patches; pch++) {
		aptr = flow_table[pch].adj_list;
		total_area += flow_table[pch].area;

		// check to see if it is a stream patch
		// if is is, add the distance to downstream patch
		if ((flow_table[pch].land == 1)) {
			n_adjacent = 0;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				inx = aptr->inx;

				if (aptr->gamma > 0.0) {
					xrun = (double)(pow((flow_table[pch].x - flow_table[inx].x), 2.0));
					yrun = (double)(pow((flow_table[pch].y - flow_table[inx].y), 2.0));
					run = (double)(sqrt(xrun + yrun) * (cell));

					drainage_density += run;
					if (flow_table[inx].land != 1)
						std::cout << "Downstream patch from is not a stream " << flow_table[pch].patchID << " " <<
						flow_table[inx].patchID << std::endl;
					n_adjacent += 1;
				}

				aptr = aptr->next;
			}

			if (n_adjacent == 0) {
				std::cout << "STREAM to outlet " << flow_table[pch].patchID << std::endl;
				drainage_density += sqrt(double(flow_table[pch].area * cell));
			}
		}
	}

	drainage_density = (double)(drainage_density / (cell * total_area));
	std::cout << "Drainage Density is " << drainage_density << std::endl;

	return;
}

//DEVIDE ALGORITHM FLAGS
void  compute_gamma(struct flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx)
{

	cout << "START COMPUTING GAMMA OF EACH PATCH AND IT'S NEIGHBORS: thread "<<thread_inx<< endl;

	if (algorithm_flag == 0)
		compute_gamma_D8(flow_table, num_patches, sc_flag, cell, algorithm_flag, pch_boun, thread_inx);
	else if (algorithm_flag == 1)
		compute_gamma_MD8(flow_table, num_patches, sc_flag, cell, algorithm_flag, pch_boun, thread_inx);
	else if (algorithm_flag == 2)
		compute_gamma_D_inf(flow_table, num_patches, sc_flag, cell, algorithm_flag, pch_boun, thread_inx);
	else if (algorithm_flag == 3 || algorithm_flag == 5)
		compute_gamma_RMD_inf(flow_table, num_patches, sc_flag, cell, algorithm_flag, pch_boun, thread_inx);
	else if (algorithm_flag == 4)
		compute_gamma_MFD_MD(flow_table, num_patches, sc_flag, cell, algorithm_flag, pch_boun, thread_inx);
	else
		cerr << "An undefined algorithm_flag." << endl;

	cout << "END OF COMPUTING GAMMAS" << endl;
}

void  compute_gamma_D8(struct flow_struct *flow_table, int num_patches, int sc_flag, double cell,int algorithm_flag, int*pch_boun,int thread_inx)
{
	// local variable declarations
	int p; //,z,h
	int inx;
	int neigh;
	int pch;
	int single_index;
	int num_acc[9]{};
	double  rise, run;
	double xrun, yrun;
	double max_slope;
	struct adj_struct *aptr;

	// compute mean pch values
	for (pch = 1; pch <= num_patches; pch++) {

		if (sc_flag == 1)
			flow_table[pch].internal_slope = flow_table[pch].internal_slope;

		flow_table[pch].flna = flow_table[pch].flna / flow_table[pch].area;
		flow_table[pch].total_gamma = 0;
		flow_table[pch].road_dist = 0;
		flow_table[pch].inflow_cnt = 0;

		flow_table[pch].slope = 0;
		flow_table[pch].real_adjacent = 0;
	}

	// calculate gamma for each neighbour
	for (pch = pch_boun[thread_inx]+1; pch <= pch_boun[thread_inx+1]; pch++) {

		aptr = flow_table[pch].adj_list;

		max_slope = 0;
		single_index = 0;

		//1. COMPUTE THE total_gamma based on Quinn, 1991
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {

			inx = aptr->inx;
			// rise is difference in elevation between two cells
			// run is distance between centorid of two neighbors
			rise = flow_table[pch].z - flow_table[inx].z;
			xrun = pow((flow_table[pch].x - flow_table[inx].x), 2.0);
			yrun = pow((flow_table[pch].y - flow_table[inx].y), 2.0);
			run = sqrt(xrun + yrun); // * (cell);
			aptr->slope = (double)(rise / run);


			aptr->gamma = (double)(aptr->perimeter *1.* aptr->slope);

			if (aptr->slope > max_slope)
				max_slope = aptr->slope;

			aptr->z = flow_table[inx].z;


			//if neighbor's elevation is higher, then aptr->gamma is less than zero, which
			//means water can not flow upslope
			if (aptr->gamma < 0.0) aptr->gamma = 0;

			//sum up all gammas
			flow_table[pch].total_gamma += aptr->gamma;

			if (aptr->slope > 0)
				flow_table[pch].slope += aptr->slope;

			aptr = aptr->next;
		}

		//2. DISTRIBUTE GAMMA TO ONLY ONE NEIGHBOUR
		aptr = flow_table[pch].adj_list;
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			if (single_index == 0 && fabs(aptr->slope - max_slope) < 0.00000000001 && pch != num_patches) {
				aptr->gamma = 1; single_index++;
			}
			else
				aptr->gamma = 0;
			aptr = aptr->next;
		}


		//3. USING D8 IN STREAM STYPE CELLS
		double min_z_stream{ 99999 };//as an infinite z value
		if (flow_table[pch].land == 1 && pch != num_patches) {

			//find min_z_stream
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (aptr->z < flow_table[pch].z && flow_table[aptr->inx].land == 1 && aptr->z < min_z_stream)
					min_z_stream = aptr->z;
				aptr = aptr->next;
			}

			//check if there is no stream neighbour
			if (fabs(min_z_stream - 99999) < 0.000000001)
			{
				cout << "there is no stream neighbour!" << endl;
				getchar();
			}

			//redivided
			single_index = 0;
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (flow_table[(aptr->inx)].land == 1) {
					if (fabs(aptr->z - min_z_stream) < 0.0000000001 && single_index == 0) {
						aptr->gamma = 1; single_index++;
					}
					else aptr->gamma = 0;
				}
				else
					aptr->gamma = 0;
				aptr = aptr->next;
			}
		}//end of rectified in stream


		// divided by total_gamma
		aptr = flow_table[pch].adj_list;
		if (flow_table[pch].total_gamma == 0.0) {
			flow_table[pch].slope = 0.0;
		}
		else
			flow_table[pch].slope = flow_table[pch].slope / flow_table[pch].total_gamma;

		//4. COMPUTING THE  NUMBER OF DOWNSLOPE CELLS
		double x = 0;
		int num_patch = 0;
		aptr = flow_table[pch].adj_list;
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			if (aptr->gamma > 0) {
				x += aptr->gamma;
				num_patch++;
				flow_table[pch].number++;
			}
			aptr = aptr->next;
		}
		num_acc[num_patch]++;

		//CHECK IF TOTAL_GAMMA IS TOO SMALL
		if (flow_table[pch].total_gamma < 0.01) flow_table[pch].total_gamma = 0.01;

		if (pch % 10000 == 0) cout << "\t" << pch / 10000 << " finished" << endl;

		if (flow_table[pch].patchID == -99999 )
			cout << "check" << endl;


	}//end of a pch
	cout << "D8\nNumber of downslope cells receiving water:"<<thread_inx << endl;
	for (int i = 0; i < 9; i++)
		cout << i << ": " << num_acc[i] << endl;

}//END OF D8

void compute_gamma_MD8(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx)
{
	int p; //,z,h
	int inx;
	int neigh;
	int pch;

	double  rise, run; //mult,
	double xrun, yrun;
	struct adj_struct *aptr;
	int num_acc[9]{};

	// compute mean pch values
	for (pch = 1; pch <= num_patches; pch++) {
		flow_table[pch].x = flow_table[pch].x;
		flow_table[pch].y = flow_table[pch].y;
		flow_table[pch].z = flow_table[pch].z;
		if (sc_flag == 1)
			flow_table[pch].internal_slope = flow_table[pch].internal_slope;

		flow_table[pch].flna = flow_table[pch].flna / flow_table[pch].area;
		flow_table[pch].total_gamma = 0.;
		flow_table[pch].road_dist = 0.;
		flow_table[pch].inflow_cnt = 0;
	}

	// calculate gamma for each neighbour
	for (pch = pch_boun[thread_inx] + 1; pch <= pch_boun[thread_inx + 1]; pch++) {
		aptr = flow_table[pch].adj_list;
		flow_table[pch].slope = 0.0;

		//1. COMPUTE THE total_gamma based on Quinn, 1991
		//HERE WE DIDN'T APPLY THE FLOW WIDTH 0.352 0.5
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {

			inx = aptr->inx;
			// rise is difference in elevation between two cells
			// run is distance between centorid of two neighbors
			rise = flow_table[pch].z - flow_table[inx].z;
			xrun = pow((flow_table[pch].x - flow_table[inx].x), 2.0);
			yrun = pow((flow_table[pch].y - flow_table[inx].y), 2.0);
			run = sqrt(xrun + yrun);
			aptr->slope = rise / run;
			aptr->gamma = aptr->perimeter * aptr->slope;

			aptr->z = flow_table[inx].z;

			if (aptr->gamma < 0.0) aptr->gamma = 0.0;

			flow_table[pch].total_gamma += aptr->gamma;
			if (aptr->slope > 0)
				flow_table[pch].slope += aptr->slope;

			aptr = aptr->next;
		}


		aptr = flow_table[pch].adj_list;
		if (flow_table[pch].total_gamma == 0.0)
			flow_table[pch].slope = 0.0;

		//2. DISTRIBUTE GAMMA TO ONLY ONE NEIGHBOUR
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			//xu. THE FLOW WIDTH SHOULD BE APPLIED IN MD8 
			if (flow_table[pch].total_gamma > 0)
				aptr->gamma = aptr->slope *aptr->perimeter / flow_table[pch].total_gamma;
			else
				aptr->gamma = 0.0;
			aptr = aptr->next;
		}

		//3. USING D8 IN STREAM STYPE CELLS 
		double min_z_stream{ 99999 };
		if (flow_table[pch].land == 1 && pch != num_patches) {
			//find min_z_stream
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (aptr->z < flow_table[pch].z && flow_table[aptr->inx].land == 1 && aptr->z < min_z_stream)
					min_z_stream = aptr->z;
				aptr = aptr->next;
			}
			//check if there is no stream neighbour
			if (fabs(min_z_stream - 99999) < 0.00001)
			{
				cout << "there is no stream neighbour!" << endl;
				getchar();
			}
			//redivided
			int num(0);
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (flow_table[aptr->inx].land == 1) {
					if (fabs(aptr->z - min_z_stream) < 0.0000000001&&num == 0) {
						aptr->gamma = 1; num++;
					}
					else aptr->gamma = 0;
				}
				else
					aptr->gamma = 0;
				aptr = aptr->next;
			}
		}//end of rectified in stream


		 //4. COMPUTING THE  NUMBER OF DOWNSLOPE CELLS
		double x = 0;
		int num_patch = 0;
		aptr = flow_table[pch].adj_list;
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			if (aptr->gamma > 0.00001) {
				x += aptr->gamma;
				num_patch++;
				flow_table[pch].number++;
			}
			aptr = aptr->next;
		}
		num_acc[num_patch]++;

		//CHECK IF TOTAL_GAMMA IS TOO SMALL
		if (flow_table[pch].total_gamma < 0.001) flow_table[pch].total_gamma = 0.001;
		if (pch % 10000 == 0) cout <<"\t"<< pch/ 10000 << " finished" << endl;
	}//end of a pch
	cout << "MD8\nNumber of downslope cells receiving water:" <<thread_inx<< endl;
	for (int i = 0; i < 9; i++)
		cout << i << ": " << num_acc[i] << endl;
	cout << "End of calcuating gamma.\n" << endl;
	return;
}

void compute_gamma_MFD_MD(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx)
{
	//rectify for sub-basin
	int rectify_flow_routing_flag = 1;

	int p; //,z,h
	int inx;
	int neigh;
	int pch;

	double  rise, run; //mult,
	double xrun, yrun;
	struct adj_struct *aptr;
	int num_acc[9]{};

	// compute mean pch values
	for (pch = 1; pch <= num_patches; pch++) {
		flow_table[pch].x = flow_table[pch].x;
		flow_table[pch].y = flow_table[pch].y;
		flow_table[pch].z = flow_table[pch].z;
		if (sc_flag == 1)
			flow_table[pch].internal_slope = flow_table[pch].internal_slope;

		flow_table[pch].flna = flow_table[pch].flna / flow_table[pch].area;
		flow_table[pch].total_gamma = 0.;
		flow_table[pch].road_dist = 0.;
		flow_table[pch].inflow_cnt = 0;
	}


	// calculate gamma for each neighbour
	for (pch = pch_boun[thread_inx] + 1; pch <= pch_boun[thread_inx + 1]; pch++) {
		aptr = flow_table[pch].adj_list;
		flow_table[pch].slope = 0.0;

		double total_slope{};
		double max_slope{};

		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
			
			inx = aptr->inx;
			// rise is difference in elevation between two cells
			// run is distance between centorid of two neighbors
			rise = flow_table[pch].z - flow_table[inx].z;
			xrun = pow((flow_table[pch].x - flow_table[inx].x), 2.0);
			yrun = pow((flow_table[pch].y - flow_table[inx].y), 2.0);
			run = sqrt(xrun + yrun); // * (cell);
			aptr->slope = rise / run;

			//FIND MAX SLOPE FOR MFD-MD
			if (aptr->slope > 1)
				max_slope = 1;
			else if (aptr->slope > max_slope)
				max_slope = aptr->slope;

			aptr->gamma = aptr->slope*aptr->perimeter;

			aptr->z = flow_table[inx].z;

			if (aptr->gamma < 0.0) aptr->gamma = 0.0;
			if (aptr->slope < 0.0) aptr->slope = 0.0;

			flow_table[pch].total_gamma += aptr->gamma;
			flow_table[pch].slope += aptr->slope;

			aptr = aptr->next;
		}

		//2. ADD EXPONENTTIAL VALUE max_slope*8.9 + 1.1
		aptr = flow_table[pch].adj_list;
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
			aptr->slope = pow(aptr->slope, max_slope*8.9 + 1.1);//A KEY EQUATION TO DIFFER MD8 AND MFD-MD
			
			if (rectify_flow_routing_flag == 1)
				if (aptr->subbasin_id != flow_table[pch].subbasin_id)
					aptr->slope = 0;

			total_slope += aptr->slope*aptr->perimeter;
			aptr = aptr->next;
		}

		// divided by total_gamma
		aptr = flow_table[pch].adj_list;

		if (flow_table[pch].total_gamma == 0.0)
			flow_table[pch].slope = 0.0;

		//3. DISTRIBUTE GAMMA TO  NEIGHBOURS BY SLOPE
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			//xU. THE FLOW WIDTH SHOULD BE APPLIED IN MFD_MD
			if (flow_table[pch].total_gamma != 0.0)
				aptr->gamma = aptr->slope*aptr->perimeter / total_slope;
			else
				aptr->gamma = 0.0;
			aptr = aptr->next;
		}

		//3. USING D8 IN STREAM STYPE CELLS 
		double min_z_stream{ 99999 };
		if (flow_table[pch].land == 1 && pch != num_patches) {
			//find min_z_stream
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (aptr->z < flow_table[pch].z && flow_table[aptr->inx].land == 1 && aptr->z < min_z_stream)
					min_z_stream = aptr->z;
				aptr = aptr->next;
			}
			//check if there is no stream neighbour
			if (fabs(min_z_stream - 99999) < 0.00001)
			{
				cout << "there is no stream neighbour!" << endl;
				getchar();
			}
			//redivided
			int num(0);
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (flow_table[aptr->inx].land == 1) {
					if (fabs(aptr->z - min_z_stream) < 0.0000000001&&num == 0) {
						aptr->gamma = 1; num++;
					}
					else aptr->gamma = 0;
				}
				else
					aptr->gamma = 0;
				aptr = aptr->next;
			}
		}//end of rectified in stream

		//4. COMPUTING THE  NUMBER OF DOWNSLOPE CELLS
		double x = 0;
		int num_patch = 0;
		aptr = flow_table[pch].adj_list;
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			if (aptr->gamma > 0.00001) {
				x += aptr->gamma;
				num_patch++;
				flow_table[pch].number++;
			}
			aptr = aptr->next;
		}
		num_acc[num_patch]++;

		//CHECK IF GAMMA IS TOO SMALL
		if (flow_table[pch].total_gamma < 0.001) flow_table[pch].total_gamma = 0.001;
		if (pch % 10000 == 0) cout << "\t" << pch / 10000 << " finished" << endl;

	}//end of a pch
	cout << "MFD-md\nNumber of downslope cells receiving water:" <<thread_inx<< endl;
	for (int i = 0; i < 9; i++)
		cout << i << ": " << num_acc[i] << endl;
	cout << "End of calcuating gamma.\n" << endl;
	return;
}

void compute_gamma_D_inf(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx)
{


	// local variable declarations
	int p; //,z,h
	int inx;
	int neigh;
	int pch;
	int total{ 0 };
	int num_acc[9]{};

	int real_adjacent;//used when neighbors is more than 2
	int inx_max1, inx_max2_a = 0, inx_max2_b = 0, inx_max2;//as the inx of neigbor with maxest slope and maxest near it
	double max_slope1 = 0, max_slope2 = 0;

	double  rise, run;
	double xrun, yrun;
	struct adj_struct *aptr;//pointer for reading neighbors

	double d1, d2;
	double g1 = 0., g2 = 0.;

	// compute mean pch values
	for (pch = 1; pch <= num_patches; pch++) {
		flow_table[pch].x = flow_table[pch].x;
		flow_table[pch].y = flow_table[pch].y;
		flow_table[pch].z = flow_table[pch].z;//just unused
		if (sc_flag == 1)//sc_flag=1 here
			flow_table[pch].internal_slope = flow_table[pch].internal_slope;

		flow_table[pch].flna = flow_table[pch].flna / flow_table[pch].area;
		flow_table[pch].total_gamma = 0.;
		flow_table[pch].road_dist = 0.;
		flow_table[pch].inflow_cnt = 0;
	}

	//Calculate gamma for each neighbour
	for (pch = pch_boun[thread_inx] + 1; pch <= pch_boun[thread_inx + 1]; pch++) 
	{
		aptr = flow_table[pch].adj_list;//initialized, after it elements of list will not be changed
		flow_table[pch].slope = 0.0;
		real_adjacent = 0;
		inx_max1 = 0;
		inx_max2_a = 0; inx_max2_b = 0; inx_max2 = 0;
		max_slope1 = 0;


		//1. COMPUTE THE total_gamma based on Quinn, 1991
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)//find the labor with maxest slope
		{
			inx = aptr->inx;
			// rise is difference in elevation between two cells
			// run is distance between centorid of two neighbors
			//here we assume there is only one neighbor
			rise = flow_table[pch].z - flow_table[inx].z;
			xrun = pow((flow_table[pch].x - flow_table[inx].x), 2.0);
			yrun = pow((flow_table[pch].y - flow_table[inx].y), 2.0);
			run = sqrt(xrun + yrun); // * (cell);

			aptr->slope = (double)(rise / run);

			aptr->gamma = (double)(aptr->perimeter *1.* aptr->slope);

			//cout << pch << " " << aptr->perimeter <<"  "<<aptr->slope<< endl;

			aptr->z = flow_table[inx].z;

			if (aptr->slope >= max_slope1)
			{
				max_slope1 = aptr->slope;
				inx_max1 = aptr->inx;
			}
			aptr->z = flow_table[inx].z;

			if (rise > 0)
				real_adjacent += 1;

			if (aptr->gamma < 0.0) aptr->gamma = 0.0;
			flow_table[pch].total_gamma += aptr->gamma;
			aptr->gamma = 0.0;

			flow_table[pch].slope += aptr->slope;
			aptr = aptr->next;
		}

		//2.COMPUTE GAMMA OF D_INF, 
		//THOUGH USING DIFFERENT ALGORITHMS TO RMD_INF MD_INF 
		//IT IS EVALUATED TO HAVE THE SAME RESULTS WITHIN EACH FACETS

		if (real_adjacent >= 2)//when >=2
		{
			aptr = flow_table[pch].adj_list;//reinitialized

			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)//find the max2
			{
				if (fabs(pow(flow_table[inx_max1].x - flow_table[aptr->inx].x, 2) + pow(flow_table[inx_max1].y - flow_table[aptr->inx].y, 2) - pow(cell, 2)) < 0.00000001 && inx_max2_a == 0)
				{
					inx_max2_a = aptr->inx; aptr = aptr->next; continue;

				}
				if (fabs(pow(flow_table[inx_max1].x - flow_table[aptr->inx].x, 2) + pow(flow_table[inx_max1].y - flow_table[aptr->inx].y, 2) - pow(cell, 2)) < 0.00000001  && inx_max2_b == 0)
				{
					inx_max2_b = aptr->inx; aptr = aptr->next; continue;
				}
				aptr = aptr->next;
			}

			aptr = flow_table[pch].adj_list;//reinitialized

			if (inx_max2_a == 0)//when there's no neighbors near the steepest path
			{
				inx_max2 = 0;
				aptr = flow_table[pch].adj_list;//reinitialized
				for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
				{
					if (aptr->inx == inx_max1)
						aptr->gamma = 1.0;
					else
						aptr->gamma = 0.;
					aptr = aptr->next;
				}

			}
			else if (inx_max2_b == 0)//when there's only neighbors near the steepest path
			{
				inx_max2 = inx_max2_a;
				aptr = flow_table[pch].adj_list;//reinitialized
				//calculating gamma for each patch
				if (fabs((flow_table[inx_max1].x - flow_table[pch].x)*(flow_table[inx_max1].y - flow_table[pch].y)) > 0.00000001)//diagonal
				{
					double sita;
					d2 = flow_table[pch].z - flow_table[inx_max2].z;
					d1 = flow_table[inx_max2].z - flow_table[inx_max1].z;
					if (d1 <= 0)
					{
						g1 = 0.;
						g2 = 1.0;
					}
					else if (d1 >= d2)
					{
						g1 = 1.0;
						g2 = 0.;
					}
					else
					{
						sita = atan(d1*1.0 / d2);
						g1 = sita * 1.0 / (0.25*PI);
						g2 = 1.0 - g1;
					}
					aptr = flow_table[pch].adj_list;//reinitialized
					for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)//distributing gamma
					{
						if (static_cast<int>(aptr->inx) == inx_max1)
							aptr->gamma = g1;
						else if (static_cast<int>(aptr->inx) == inx_max2)
							aptr->gamma = g2;
						else
							aptr->gamma = 0.;
						aptr = aptr->next;
					}
					aptr = flow_table[pch].adj_list;//reinitialized
				}

				else//vertical or horizontal
				{
					double d1, d2;//d1 for d between pch.z-max1.z, d2 for d between max1.z -max2.z
					double g1, g2;
					double sita = 0;
					d2 = flow_table[pch].z - flow_table[inx_max1].z;
					d1 = flow_table[inx_max1].z - flow_table[inx_max2].z;
					if (d1 <= 0)
					{
						g1 = 1.0;
						g2 = 0.;
					}
					else if (d1 >= d2)
					{
						g1 = 0.;
						g2 = 1.0;
					}
					else
					{
						sita = atan(d1*1.0 / d2);
						g2 = sita * 1.0 / (0.25*PI);
						g1 = 1.0 - g2;
					}
					aptr = flow_table[pch].adj_list;//reinitialized
					for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)//distributing gamma
					{
						if (aptr->inx == inx_max1)
							aptr->gamma = g1;
						else if (aptr->inx == inx_max2)
						{
							aptr->gamma = g2;
						}
						else
							aptr->gamma = 0.;
						aptr = aptr->next;
					}
				}
			}

			else //when there's two neighbors near the steepest path
			{
				if (flow_table[inx_max2_a].z >= flow_table[inx_max2_b].z)
					inx_max2 = inx_max2_b;
				else
					inx_max2 = inx_max2_a;
				aptr = flow_table[pch].adj_list;//reinitialized
				//calculating gamma for each patch
				if (fabs((flow_table[inx_max1].x - flow_table[pch].x)*(flow_table[inx_max1].y - flow_table[pch].y)) > 0.00000001)//diagonal
				{
					double sita;
					d2 = flow_table[pch].z - flow_table[inx_max2].z;
					d1 = flow_table[inx_max2].z - flow_table[inx_max1].z;

					if (d1 <= 0)
					{
						g1 = 0.;
						g2 = 1.0;
					}
					else if (d1 >= d2)
					{
						g1 = 1.0;
						g2 = 0.;
					}
					else
					{
						sita = atan(d1*1.0 / d2);
						g1 = sita * 1.0 / (0.25*PI);
						g2 = 1.0 - g1;
					}
					aptr = flow_table[pch].adj_list;//reinitialized
					for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)//distributing gamma
					{
						if (static_cast<int>(aptr->inx) == inx_max1)
							aptr->gamma = g1;
						else if (static_cast<int>(aptr->inx) == inx_max2)
							aptr->gamma = g2;
						else
							aptr->gamma = 0.;
						aptr = aptr->next;
					}
					aptr = flow_table[pch].adj_list;//reinitialized
				}

				else//vertical or horizontal
				{
					double d1, d2;//d1 for d between pch.z-max1.z, d2 for d between max1.z -max2.z
					double g1, g2;
					double sita = 0;
					d2 = flow_table[pch].z - flow_table[inx_max1].z;
					d1 = flow_table[inx_max1].z - flow_table[inx_max2].z;
					if (d1 <= 0)
					{
						g1 = 1.0;
						g2 = 0.;
					}
					else if (d1 >= d2)
					{
						g1 = 0.;
						g2 = 1.0;
					}
					else
					{
						sita = atan(d1*1.0 / d2);
						g2 = sita * 1.0 / (0.25*PI);
						g1 = 1.0 - g2;
					}
					aptr = flow_table[pch].adj_list;//reinitialized
					for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)//distributing gamma
					{
						if (aptr->inx == inx_max1)
							aptr->gamma = g1;
						else if (aptr->inx == inx_max2)
						{
							aptr->gamma = g2;
						}
						else
							aptr->gamma = 0.;
						aptr = aptr->next;
					}
				}//end of vertical
			}//end of two neighbors near the steepest path
		}//end of two real adjacecent
		else if (real_adjacent == 1)//WHEN THERE ARE ONLY ONE DOWNSLOPE NEIGHBOURS
		{

			aptr = flow_table[pch].adj_list;//reinitialized
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
			{
				if (aptr->inx == inx_max1)
					aptr->gamma = 1.0;
				else
					aptr->gamma = 0.;
				aptr = aptr->next;
			}
		}

		//3.RECTIFIED FOR OUTLET OF A BASIN
		if (pch == num_patches)
		{
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
			{
				aptr->gamma = 0.;
				aptr = aptr->next;
			}
		}

		//3. USING D8 IN STREAM STYPE CELLS 
		double min_z_stream{ 99999 };
		if (flow_table[pch].land == 1 && pch != num_patches) {
			//find min_z_stream
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (aptr->z < flow_table[pch].z && flow_table[aptr->inx].land == 1 && aptr->z < min_z_stream)
					min_z_stream = aptr->z;
				aptr = aptr->next;
			}
			//check if there is no stream neighbour
			if (fabs(min_z_stream - 99999) < 0.00001)
			{
				cout << "there is no stream neighbour!" << endl;
				getchar();
			}
			//redivided
			int num(0);
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (flow_table[aptr->inx].land == 1) {
					if (fabs(aptr->z - min_z_stream) < 0.0000000001&&num == 0) {
						aptr->gamma = 1; num++;
					}
					else aptr->gamma = 0;
				}
				else
					aptr->gamma = 0;
				aptr = aptr->next;
			}
		}//end of rectified in stream

		 //4. COMPUTING THE  NUMBER OF DOWNSLOPE CELLS
		double x = 0;
		int num_patch = 0;
		aptr = flow_table[pch].adj_list;
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			if (aptr->gamma > 0.00001) {
				x += aptr->gamma;
				num_patch++;
				flow_table[pch].number++;
			}
			aptr = aptr->next;
		}
		num_acc[num_patch]++;
		if (pch % 10000 == 0) cout << "\t" << pch / 10000 << " finished" << endl;

	}//end of a pch
	cout << "D_inf\nNumber of downslope cells receiving water:" <<thread_inx<< endl;
	for (int i = 0; i < 9; i++)
		cout << i << ": " << num_acc[i] << endl;
	cout << "End of calcuating gamma.\n" << endl;
}

//THIS FUNCTION ARE USED IN BOTH RMD_INF AND MD_INF, CONTROLED BY Algorithm_flag
void compute_gamma_RMD_inf(flow_struct *flow_table, int num_patches, int sc_flag, double cell, int algorithm_flag, int*pch_boun, int thread_inx)
{
	// local fuction declarations
	char algorithm[6][10]{ "D8", "MD8", "D_inf", "RMD_inf", "MFD_md", "MD_inf" };
	void computing_sk(struct flow_struct *flow_table, struct adj_struct* ptr1, int pch, int acc_num, double cell);
	int MD_inf_acc(flow_struct *flow_table, int pch, int i_val_arr[8], int acc_i[8][9], int acc_num[8]);

	// local variable declarations
	int p; //,z,h
	int inx;
	int neigh;
	int pch;

	int acc_i[8][9]{};
	int acc_num[8]{};
	int routing_patch_num[9]{};


	double x_value, y_value;
	int s = 0;

	double  rise, run;
	double xrun, yrun;
	struct adj_struct *aptr;//pointer for reading neighbors
	struct adj_struct *aptr1;

	struct adj_struct *ptr1, *ptr2;//used in calculation within a tri


	// compute mean pch values
	for (pch = 1; pch <= num_patches; pch++) {
		flow_table[pch].x = flow_table[pch].x;
		flow_table[pch].y = flow_table[pch].y;
		flow_table[pch].z = flow_table[pch].z;//just unused
		if (sc_flag == 1)//sc_flag=1 here
			flow_table[pch].internal_slope = flow_table[pch].internal_slope;

		flow_table[pch].flna = flow_table[pch].flna / flow_table[pch].area;
		flow_table[pch].total_gamma = 0.;
		flow_table[pch].road_dist = 0.;
		flow_table[pch].inflow_cnt = 0;
	}

	// calculate gamma for each neighbour
	for (pch = pch_boun[thread_inx] + 1; pch <= pch_boun[thread_inx + 1]; pch++) 
	{
		//cout << "\n pch is " << pch << "\n";
		aptr = flow_table[pch].adj_list;//initialized, after it elements of list will not be changed
		flow_table[pch].slope = 0.0;
		flow_table[pch].real_adjacent = 0;


		//1. COMPUTING I_VALUE AND SOME INDEX
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)//initialized all data
		{
			inx = aptr->inx;

			//distributing i_value, for sorting acc_area in this small region
			// 1   2   3
			// 8  pch  4
			// 7   6   5
			x_value = flow_table[pch].x - flow_table[inx].x;
			y_value = flow_table[pch].y - flow_table[inx].y;
			if (fabs(x_value) < 0.0000001)
			{
				if (fabs(y_value - cell) < 0.0000001)
					aptr->i_value = 2;
				else if (fabs(y_value + cell) < 0.00000001)
					aptr->i_value = 6;
			}
			else if (fabs(x_value - cell) < 0.0000001)
			{
				if (fabs(y_value - cell) < 0.0000001)
					aptr->i_value = 1;
				else if (fabs(y_value) < 0.0000001)
					aptr->i_value = 8;
				else if (fabs(y_value + cell) < 0.0000001)
					aptr->i_value = 7;
			}
			else if (fabs(x_value + cell) < 0.0000001)
			{
				if (fabs(y_value - cell) < 0.0000001)
					aptr->i_value = 3;
				else if (fabs(y_value) < 0.0000001)
					aptr->i_value = 4;
				else if (fabs(y_value + cell) < 0.0000001)
					aptr->i_value = 5;
			}
			// rise is difference in elevation between two cells
			// run is distance between centorid of two neighbors
			//here we assume there is only one neighbor
			rise = flow_table[pch].z - flow_table[inx].z;
			xrun = pow((flow_table[pch].x - flow_table[inx].x), 2.0);
			yrun = pow((flow_table[pch].y - flow_table[inx].y), 2.0);
			run = sqrt(xrun + yrun); // * (cell);
			aptr->slope = (double)(rise / run);

			aptr->gamma = 0;//initialized as 0.

			if (aptr->slope > 0.0)
				flow_table[pch].total_gamma += aptr->perimeter * aptr->slope;

			aptr->z = flow_table[inx].z;


			flow_table[pch].slope += aptr->slope;
			aptr = aptr->next;
		}//END OF 1

		//2. SORTING I_VALUES LOWER
		aptr = flow_table[pch].adj_list;
		int i_val_arr[9]{};//i_value collection
		int t{};
		for (int t = 0; t < flow_table[pch].num_adjacent; t++)//initialized
		{
			i_val_arr[t] = aptr->i_value;
			aptr = aptr->next;
		}
		for (int i = 0; i < flow_table[pch].num_adjacent; i++)//sorting
		{
			for (int j = flow_table[pch].num_adjacent - 1; j > i; j--) {
				if (i_val_arr[j] < i_val_arr[i]) {
					int mid{};
					mid = i_val_arr[j];
					i_val_arr[j] = i_val_arr[i];
					i_val_arr[i] = mid;
				}
			}

		}//end of sorting

		//3. Distributing i_value for MD_inf and RMD_inf
		//A KEY FUNCTION FOR DEPARTING ACC
		flow_table[pch].num_acc = MD_inf_acc(flow_table, pch, i_val_arr, acc_i, acc_num);

		//Deliver ACC i_value num
		for (neigh = 0; neigh < 8; neigh++) {
			flow_table[pch].acc_patch_num[neigh] = acc_num[neigh];
		}

		//4. BUILDING LISTS FOR ALL ACC
		struct adj_struct *p[8][9];

		for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++) {
			p[neigh][0] = (struct adj_struct *)malloc(sizeof(struct adj_struct));
			flow_table[pch].adj_list_acc[neigh] = p[neigh][0];
			for (t = 0; t < 8; t++)
			{
				p[neigh][t + 1] = (struct adj_struct *)malloc(sizeof(struct adj_struct));
				p[neigh][t]->next = p[neigh][t + 1];
			}
		}

		//INITIAL LIST
		for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
		{
			aptr1 = flow_table[pch].adj_list_acc[neigh];
			for (int t = 0; t < flow_table[pch].acc_patch_num[neigh]; t++)
			{
				//FIND VALUE
				aptr = flow_table[pch].adj_list;
				for (int j = 1; j <= flow_table[pch].num_adjacent; j++)
				{
					if (acc_i[neigh][t] == aptr->i_value) break;
					aptr = aptr->next;
				}
				aptr1->patchID = aptr->patchID;
				aptr1->landID = aptr->landID;
				aptr1->inx = aptr->inx;
				aptr1->perimeter = aptr->perimeter;
				aptr1->slope = aptr->slope;
				aptr1->gamma = 0.;
				aptr1->z = aptr->z;
				aptr1->i_value = aptr->i_value;
				aptr1->DFA = 0;
				aptr1->DFB = 0;
				aptr1->DF = 0;
				aptr1->DF_ratio = 0;
				aptr1->distribute_ratio = 0;
				aptr1->max_slope = 0;
				aptr1 = aptr1->next;
			}
		}

		//5. Sk and Maximun slope
		if (pch != num_patches && flow_table[pch].num_acc != 0)
		{
			//within every tri
			for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++) {
				ptr1 = flow_table[pch].adj_list_acc[neigh];
				computing_sk(flow_table, ptr1, pch, neigh, cell);
			}
		}// end of all Sk

		//6. Tk
		//FOR TOTAL SLOPE
		double total_slope = 0.;
		flow_table[pch].DF_num = 0;
		//DESIDE DF
		for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
		{
			ptr1 = flow_table[pch].adj_list_acc[neigh];
			if (flow_table[pch].acc_patch_num[neigh] == 1) {
				if (ptr1->DF == 1) {
					total_slope += ptr1->max_slope;
					flow_table[pch].DF_num++;
				}
			}
			else {
				for (int i = 0; i < flow_table[pch].acc_patch_num[neigh] - 1; i++) {

					//AN IMPORTANT ISSUE DIFFER BETWEEN THESE TWO ALGORITHMS
					//FIX FOR RMD_inf
					if (algorithm_flag == 3) {
						if (ptr1->z > flow_table[pch].z  && ptr1->next->z > flow_table[pch].z) {
							ptr1->DF = 0;
						}
						else if (ptr1->z > flow_table[pch].z || ptr1->next->z > flow_table[pch].z) {
							ptr1->DF = 1;
							ptr1->max_slope = ptr1->max_slope*0.5;//SET IT TO BE O.5
						}
						else
							ptr1->DF = 1;
					}
					else if (algorithm_flag == 5) {
						if (flow_table[pch].acc_patch_num[neigh] == 9 && i == 7) {
							if (ptr1->next->DFA == 1 && flow_table[pch].adj_list_acc[neigh]->DFA == 1)
								ptr1->DF = 1;
						}

						if (ptr1->DFA == 1 && ptr1->DFB == 1) {
							ptr1->DF = 1;

						}
					}
					else cout << "\nFATAL ERROR IN ALGORITHM'S FLAG\n" << endl;

					if (ptr1->DF == 1) {
						total_slope += ptr1->max_slope;
						flow_table[pch].DF_num++;
					}
					ptr1 = ptr1->next;
				}
			}
		}//end of distributing all gamma

		//Fix for no locally steepest slope
		if (flow_table[pch].DF_num < 1) {

			for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
			{
				ptr1 = flow_table[pch].adj_list_acc[neigh];
				if (flow_table[pch].acc_patch_num[neigh] == 1) {
					if (ptr1->max_slope > total_slope) total_slope = ptr1->max_slope;
				}
				else {
					for (int i = 0; i < flow_table[pch].acc_patch_num[neigh] - 1; i++) {
						if (ptr1->max_slope > total_slope) total_slope = ptr1->max_slope;
						ptr1 = ptr1->next;
					}
				}
			}
		}
		if (flow_table[pch].DF_num < 1) {
			for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
			{
				ptr1 = flow_table[pch].adj_list_acc[neigh];
				if (flow_table[pch].acc_patch_num[neigh] == 1) {
					if (fabs(ptr1->max_slope - total_slope) < 0.00000000001) { ptr1->DF = 1; break; }
				}
				else {
					for (int i = 0; i < flow_table[pch].acc_patch_num[neigh] - 1; i++) {
						if (fabs(ptr1->max_slope - total_slope) < 0.00000000001) { ptr1->DF = 1; break; }
						ptr1 = ptr1->next;
					}
				}
			}
		}

		//DEVIDE DF_RATIO AS FINAL DIVITION TO A FACET
		for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
		{
			ptr1 = flow_table[pch].adj_list_acc[neigh];
			if (flow_table[pch].acc_patch_num[neigh] == 1) {
				if (ptr1->DF == 1)
					ptr1->DF_ratio = ptr1->max_slope / total_slope;
			}
			else {
				for (int i = 0; i < flow_table[pch].acc_patch_num[neigh] - 1; i++) {
					if (ptr1->DF == 1) {
						ptr1->DF_ratio = ptr1->max_slope / total_slope;
					}
					ptr1 = ptr1->next;
				}
			}
		}//end of distributing all gamma

		//7. FINAL GAMMA IN EACH FACET
		for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
		{
			ptr1 = flow_table[pch].adj_list_acc[neigh];
			if (flow_table[pch].acc_patch_num[neigh] == 1) {
				if (ptr1->DF == 1) {
					ptr1->gamma += ptr1->DF_ratio * (ptr1->distribute_ratio);
				}
			}
			else {
				ptr2 = ptr1->next;
				for (int i = 0; i < flow_table[pch].acc_patch_num[neigh] - 1; i++)
				{
					if (ptr1->DF == 1) {
						ptr1->gamma += ptr1->DF_ratio * (ptr1->distribute_ratio);
						ptr2->gamma += ptr1->DF_ratio *  (1.0 - ptr1->distribute_ratio)*1.0;
					}
					ptr1 = ptr1->next;
					ptr2 = ptr2->next;
				}
			}
		}//end of distributing all gamma

		//8. DILIVER AND SUM UP ALL GAMMA
		for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
		{
			ptr1 = flow_table[pch].adj_list_acc[neigh];
			for (int i = 0; i < flow_table[pch].acc_patch_num[neigh]; i++)
			{
				ptr2 = flow_table[pch].adj_list;
				for (t = 1; t <= flow_table[pch].num_adjacent; t++)
				{
					if (ptr1->i_value == ptr2->i_value)
					{
						ptr2->gamma += ptr1->gamma;
					}
					ptr2 = ptr2->next;
				}
				ptr1 = ptr1->next;
			}
		}//end of distributing all gamma

		//RECTIFIED FOR OUTLET
		if (pch == num_patches)
		{
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
			{
				aptr->gamma = 0;
				aptr = aptr->next;
			}
		}

		//9. CHECK IF ITS LESS THAN 1
		double total_gamma = 0.;
		int num{};
		aptr = flow_table[pch].adj_list;
		for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++)
		{
			if (aptr->gamma < -0.00000000001) aptr->gamma = 0;
			total_gamma += aptr->gamma;
			if (aptr->gamma > 0.)  num++;//TO LAST STEP TO SUM UP
			aptr = aptr->next;
		}

		//10. USING D8 IN STREAM STYPE CELLS 
		double min_z_stream{ 99999 };
		if (flow_table[pch].land == 1 && pch != num_patches) {
			//find min_z_stream
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (aptr->z < flow_table[pch].z && flow_table[aptr->inx].land == 1 && aptr->z < min_z_stream)
					min_z_stream = aptr->z;
				aptr = aptr->next;
			}
			//check if there is no stream neighbour
			if (fabs(min_z_stream - 99999) < 0.00001)
			{
				cout << "there is no stream neighbour!" << endl;
				getchar();
			}
			//redivided
			int num(0);
			aptr = flow_table[pch].adj_list;
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				if (flow_table[aptr->inx].land == 1) {
					if (fabs(aptr->z - min_z_stream) < 0.0000000001&&num == 0) {
						aptr->gamma = 1; num++;
					}
					else aptr->gamma = 0;
				}
				else
					aptr->gamma = 0;
				aptr = aptr->next;
			}
		}//end of rectified in stream

		routing_patch_num[num]++;

		//FREE MEMORY
		for (neigh = 0; neigh < flow_table[pch].num_acc; neigh++)
			for (t = 0; t < 8; t++)
			{
				free(p[neigh][t]);
			}


		if (flow_table[pch].num_acc >= 2)
			s += 1;
		if (flow_table[pch].total_gamma == 0.0)
			flow_table[pch].slope = 0.0;
		else
			flow_table[pch].slope = flow_table[pch].slope / flow_table[pch].total_gamma;

		if (pch % 10000 == 0) cout << "\t" << pch / 10000 << " finished" << endl;

	}
	cout << algorithm[algorithm_flag] << "\nNumber of downslope cells receiving water:" <<thread_inx<< endl;
	for (int i = 0; i < 9; i++) {
		cout << "ROUTING PATCH NUM " << i << "\t" << routing_patch_num[i] << endl;
	}

	cout << "End of calculating gamma.\n";
	getchar();
}

void computing_sk(struct flow_struct *flow_table, struct adj_struct* aptr, int pch, int acc_num, double cell) {
	struct adj_struct  *ptr1, *ptr2;
	int neigh{};
	double sita{}, d1{}, d2{}, g1{}, g2{};//g1 stands for diagonal direction, d2 stands for cardinal width

	ptr1 = aptr;
	ptr2 = ptr1->next;
	flow_table[pch].acc_max_slope[acc_num] = 0.;

	if (flow_table[pch].acc_patch_num[acc_num] == 1)// while there is only one neighbor cell
	{
		if (ptr1->z < flow_table[pch].z) {
			ptr1->DF = 1;
		}
		else {
			ptr1->DF = 0;
		}
		ptr1->distribute_ratio = 1;
		ptr1->max_slope = ptr1->slope;
	}
	else {

		for (neigh = 1; neigh < flow_table[pch].acc_patch_num[acc_num]; neigh++)
		{
			if (ptr1->z < flow_table[pch].z || ptr2->z < flow_table[pch].z) {
				if (ptr1->i_value % 2 == 1)//diagonal 1 3 5 7
				{
					d2 = flow_table[pch].z - ptr2->z;
					d1 = ptr2->z - ptr1->z;
					if (d1 < 0)
					{
						g1 = 0;
						g2 = 1;
						ptr1->max_slope = d2 / cell;
						if (ptr1->max_slope < 0) cout << pch << endl;
						//DEVIDED INTO A, B
						if (ptr2->DFA == 1)
							ptr2->DFB = 1;
						else
							ptr2->DFA = 1;
					}
					else if (d1 > d2)
					{
						g1 = 1;
						g2 = 0;
						ptr1->max_slope = (d1 + d2) / (sqrt(2)*cell);
						if (ptr1->max_slope < 0) cout << pch << endl;
						//DEVIDED INTO A, B
						if (ptr1->DFA == 1)
							ptr1->DFB = 1;
						else
							ptr1->DFA = 1;

					}
					else
					{
						sita = atan(d1*1.0 / d2);
						g1 = sita * 1.0 / (0.25 * PI);//a circle is 2PI, AND 1/8
						g2 = 1.0 - g1;
						ptr1->max_slope = sqrt(d1*d1 + d2 * d2) / cell;
						ptr1->DF = 1;
					}
					ptr1->distribute_ratio = g1;
				}
				else//cardinal 2 4 6 8
				{
					d2 = flow_table[pch].z - ptr1->z;
					d1 = ptr1->z - ptr2->z;
					if (d1 < 0)
					{
						g1 = 0;
						g2 = 1;
						ptr1->max_slope = d2 / cell;

						//DEVIDED INTO A, B
						if (ptr1->DFA == 1)
							ptr1->DFB = 1;
						else
							ptr1->DFA = 1;

					}
					else if (d1 > d2)
					{
						g1 = 1;
						g2 = 0;
						ptr1->max_slope = (d1 + d2) / (sqrt(2)*cell);
						//DEVIDED INTO A, B

						if (ptr2->DFA == 1)
							ptr2->DFB = 1;
						else
							ptr2->DFA = 1;
					}
					else
					{
						sita = atan(d1*1.0 / d2);
						g1 = sita * 1.0 / (0.25 * PI);
						g2 = 1.0 - g1;
						ptr1->max_slope = sqrt(d1*d1 + d2 * d2) / cell;
						ptr1->DF = 1;
					}
					ptr1->distribute_ratio = g2;
				}
			}
			else
				ptr1->max_slope = 0;
			//start next circle

			if (ptr1->max_slope < 0.00001) {
				ptr1->max_slope = 0; ptr1->DF = 0;
			}

			ptr1 = ptr1->next;
			ptr2 = ptr2->next;

		}//finished first one
	}//end of more than 1 in triangular facet

}//end of computing sk 

//THE FIRST 
int MD_inf_acc(flow_struct *flow_table, int pch, int i_val_arr[8], int acc_i[8][9], int acc_num[8]) {

	//Extern the value
	int neigh{}, acc{}, i{}, j{};
	int acc_flag{};

	//Initial as 0
	for (acc = 0; acc < 8; acc++) {
		for (i = 0; i < 9; i++) {
			acc_i[acc][i] = 0;
		}
		acc_num[acc] = 0;
	}


	//Reinitial acc and to start distribution
	acc = 0; i = 0;
	for (neigh = 0; neigh < 8; neigh++) {

		//Stop when it's 0
		if (i_val_arr[neigh] == 0) break;

		//initial the first data
		if (neigh == 0) { acc_i[acc][i] = i_val_arr[neigh]; i++; continue; }


		//distribute i
		if (i_val_arr[neigh] != acc_i[acc][i - 1] + 1) {
			acc += 1; i = 0;
			acc_i[acc][i] = i_val_arr[neigh]; i += 1;
		}
		else {
			acc_i[acc][i] = i_val_arr[neigh];
			i += 1;
		}

		//when there is 8 and 1
		if (i_val_arr[neigh] == 8 && acc_i[0][0] == 1) {
			if (acc_i[0][7] == 8) acc_i[0][8] = 1;
			else {
				int a[8]{};
				//preserve value 
				for (int j = 0; j < 8; j++) {
					if (acc_i[0][j] != 0) {
						a[j] = acc_i[0][j];
					}
				}

				//for new array
				int ii{};
				for (int j = 0; j < 8; j++) {
					if (acc_i[acc][j] != 0) { acc_i[0][j] = acc_i[acc][j]; acc_i[acc][j] = 0; }
					else if (a[ii] != 0) { acc_i[0][j] = a[ii]; ii++; }
					else break;
				}
				acc--;//should not use acc anymore
				break;
			}
		}
	}
	//END OF SORTING 


	//3. FINAL CHECK NUM
	acc = 0;
	for (int i = 0; i < 8; i++) {
		acc_num[i] = 0;
		int acc_flag = 0;
		for (int j = 0; j < 9; j++)
		{
			if (acc_i[i][j] > 0)
			{
				if (acc_flag == 0) { acc++; acc_flag = 1; }
				acc_num[i]++;
			}

		}
	}

	//HERE WE COUT TO CHECK DATA
	/*
	if (acc > 0){
		cout << "HERE IS \t" << pch << endl;
		for (int i = 0; i < 8; i++)
		for (int j = 0; j < 9; j++)
		{
			cout << acc_i[i][j] << "\t";
			if (j == 8) cout << endl;
		}

		cout << acc << endl;
		for (int j = 0; j < acc; j++)
			cout << acc_num[j] << endl;
		//getchar();
	}
	*/

	return acc;
}

void CreatFlowTable::compute_upslope_area(flow_struct *flow_table, int num_patches, int rflag, double cell)
{
	// local fuction declarations
	//struct ID_struct  sort_flow_table();
	//int	find_patch();

	// local variable declarations
	int inx;
	int neigh;
	int pch;

	adj_struct *aptr;


	for (pch = 1; pch <= num_patches; pch++) {
		aptr = flow_table[pch].adj_list;
		if (flow_table[pch].area > 0)
			flow_table[pch].acc_area += flow_table[pch].area;
		//if(flow_table[pch].patchID==15157) {
		//	printf_s("\n %d %f \n",flow_table[pch].patchID,flow_table[pch].acc_area);
		//	getchar();
		//}
		// process roads as all accumulated area going to stream (or receiving patch)
		if ((flow_table[pch].land == 2) && (rflag == 1)) {

			inx = flow_table[pch].stream_inx;
			if (flow_table[pch].acc_area > 0)
				flow_table[inx].acc_area += flow_table[pch].acc_area;
		}

		//otherwise distribute accumulated area to downstream neighbours
		else {
			for (neigh = 1; neigh <= flow_table[pch].num_adjacent; neigh++) {
				inx = aptr->inx;
				if (aptr->gamma > 0) {
					flow_table[inx].acc_area += (flow_table[pch].acc_area) * aptr->gamma;
				}
				else
					aptr->gamma = 0;
				aptr = aptr->next;
			}
		}
	}
	return;
}

void CreatFlowTable::find_max_flna(flow_struct *flow_table, int curr, int *str_inx)
{
	adj_struct *aptr;
	int 	i, inx;
	double  max_flna;

	aptr = flow_table[curr].adj_list;
	max_flna = flow_table[aptr->inx].flna;

	// check to see if we are at the edge of the pit, find min elevation of
	// apixel draining from the pit 										

	i = 1;
	while ((i <= flow_table[curr].num_adjacent)) {
		inx = aptr->inx;
		if ((flow_table[inx].flna >= max_flna) && (aptr->gamma != 0.0)) {
			*str_inx = inx;
			max_flna = flow_table[inx].flna;
		}
		aptr = aptr->next;
		i += 1;
	}

	return;
}

void CreatFlowTable::find_min_flna(flow_struct *flow_table, int curr, int *str_inx)
{
	adj_struct *aptr;
	int 	i, inx;
	double  min_flna;

	aptr = flow_table[curr].adj_list;
	min_flna = flow_table[aptr->inx].flna;

	// check to see if we are at the edge of the pit, find min elevation of 
	// apixel draining from the pit 										

	i = 1;
	while ((i <= flow_table[curr].num_adjacent)) {
		inx = aptr->inx;

		if ((flow_table[inx].flna <= min_flna) && (aptr->gamma != 0.0)) {
			*str_inx = inx;
			min_flna = flow_table[inx].flna;
		}
		aptr = aptr->next;
		i += 1;
	}
	return;
}

int	find_patch(int num_patches, flow_struct *flow_table, int patchID)
{
	int fnd{}, inx{};

	while ((fnd == 0) && (inx <= num_patches)) {
		if (flow_table[inx].patchID == patchID)
			fnd = inx;
		else
			inx += 1;
	}

	return(fnd);
}

int CreatFlowTable::find_stream(flow_struct *flow_table, int curr, int *str_inx)
{
	adj_struct *aptr;
	long int i, inx, fnd;

	fnd = 0;
	aptr = flow_table[curr].adj_list;

	// check to see if we are at the edge of the pit, find min elevation of
	// apixel draining from the pit 									
	if (flow_table[curr].land != 1) {
		i = 1;
		while ((i <= flow_table[curr].num_adjacent) && (fnd == 0)) {
			inx = aptr->inx;

			if (aptr->gamma > 0.0)
				fnd = (int)find_stream(flow_table, inx, str_inx);

			aptr = aptr->next;
			i += 1;
		}
	}
	else {
		fnd = 1;
		*str_inx = curr;
	}

	return(fnd);
}

double	CreatFlowTable::find_top(flow_struct *flow_table, int curr, double pit_elev, int *num_pit, int *upslope_list, int *edge_inx)
{
	//int in_list(int, int *, int);
	adj_struct *aptr;
	int 	i, inx, nnew;
	double	top_elev, min_elev;
	int	next_edge_inx;

	min_elev = 0.0;
	top_elev = 0.0;

	aptr = flow_table[curr].adj_list;
	// check to see if we are at the edge of the pit, find min elevation of
	// apixel draining from the pit 										

	if ((*edge_inx == 0.0) || (flow_table[curr].z < flow_table[*edge_inx].z)) {
		for (i = 1; i <= flow_table[curr].num_adjacent; i++) {
			inx = aptr->inx;

			nnew = in_list(inx, upslope_list, *num_pit);
			if ((aptr->gamma > 0) && (nnew == 0)) {
				top_elev = flow_table[curr].z;
				if (((top_elev < min_elev) || (min_elev == 0.0)) && (top_elev > 0.0) && (aptr->z < pit_elev)) {
					*edge_inx = aptr->inx;
					min_elev = top_elev;
				}
			}

			else if ((aptr->gamma == 0) && (nnew == 0)) {
				*num_pit += 1;
				upslope_list[*num_pit] = inx;
				next_edge_inx = *edge_inx;
				top_elev = find_top(flow_table, inx, pit_elev, num_pit, upslope_list, &(next_edge_inx));
				if (((top_elev < min_elev) || (min_elev == 0.0)) && (top_elev > 0.0)) {
					*edge_inx = next_edge_inx;
					min_elev = top_elev;
				}
			}
			aptr = aptr->next;
		} // end first pass
	} //

	return (min_elev);
}

int	CreatFlowTable::in_list(int inx, int *list, int num_pit)
{
	int i, fnd;

	fnd = 0;
	for (i = 1; i <= num_pit; i++) {
		if (inx == list[i])
			fnd = i;
	}

	return(fnd);
}

double	CreatFlowTable::path_lengths(flow_struct *flow_table, int curr, int *str_inx)
{
	adj_struct  *aptr;
	int 	    i, inx, fnd;
	double      dist_neighbour;

	fnd = 0;
	aptr = flow_table[curr].adj_list;

	// check to see if we are at the edge of the pit, find min elevation of
	// apixel draining from the pit 										
	if (flow_table[curr].land != 1) {
		std::cout << "\n Current patch is %d" << flow_table[curr].patchID << std::endl;
		std::cout << "\n Next patch is %d" << flow_table[aptr->inx].patchID << std::endl;

		i = 1;
		while (i <= flow_table[curr].num_adjacent) {
			inx = aptr->inx;

			dist_neighbour = pow(30 * (flow_table[curr].x - flow_table[inx].x), 2.0) +
				pow(30 * (flow_table[curr].y - flow_table[inx].y), 2.0) +
				pow((flow_table[curr].z - flow_table[inx].z), 2.0);

			dist_neighbour = pow(dist_neighbour, 0.5);

			flow_table[inx].path_length += aptr->gamma * dist_neighbour;

			if ((flow_table[inx].path_length == 0.0) && (flow_table[inx].land != 1))
				flow_table[curr].path_length += aptr->gamma * path_lengths(flow_table, inx, str_inx);
			else
				flow_table[curr].path_length += aptr->gamma * flow_table[inx].path_length;

			aptr = aptr->next;
			i += 1;
		}
		flow_table[curr].path_length = flow_table[curr].path_length / flow_table[curr].total_gamma;
	} /* end if */
	else
		flow_table[curr].path_length = 0.0;

	return(flow_table[curr].path_length);
}

void  CreatFlowTable::remove_pits(flow_struct *flow_table, int num_patches, int sc_flag, double cell, std::ofstream & f1)
{
	// local fuction declarations
	//struct ID_struct  sort_flow_table();
	//int	   find_patch();
	//double find_top(flow_struct *, int, double, int *, int *, int *);
	//void   adjust_pit(flow_struct *, int, int, double, double);

	// local variable declarations
	int pch;
	int num_pit, num_in_pit;
	int	edge_inx;
	int *upslope_list;
	double edge_elev, top_elev;

	// allocate list
	upslope_list = new int[num_patches];
	num_pit = 0;

	for (pch = 1; pch <= num_patches; pch++) {
		// check to see if it is a pit
		if ((flow_table[pch].total_gamma == 0.) && ((flow_table[pch].land != 1))) {
			f1 << flow_table[pch].patchID << std::endl;
			num_pit += 1;
			num_in_pit = 1;
			upslope_list[1] = pch;
			edge_inx = 0;

			top_elev = find_top(flow_table, pch, flow_table[pch].z,
				&num_in_pit, upslope_list, &edge_inx);

			edge_elev = flow_table[edge_inx].z;

			if (top_elev != 0.0) {
				adjust_pit(flow_table, pch, edge_inx, edge_elev, cell);
			}
		}
	}
	//std::cout << "\n Number of pits " << num_pit << std::endl;
	f1 << " Number of pits are " << num_pit << std::endl;
	f1.close();

	//free(upslope_list)
	return;
}

void CreatFlowTable::route_roads_to_patches(flow_struct *flow_table, int num_patches, double cell, int fl_flag)
{
	//void find_min_flna(flow_struct *, int, int *);
	//int	 find_mmax_flna();

	/* local variable declarations */
	int        inx;
	int        pch;
	adj_struct *aptr;

	// calculate gamma for each neighbour
	for (pch = 1; pch <= num_patches; pch++) {
		aptr = flow_table[pch].adj_list;

		if (flow_table[pch].land == 2) {
			if (fl_flag)
				find_min_flna(flow_table, pch, &inx);
			else
				find_max_flna(flow_table, pch, &inx);

			flow_table[pch].stream_ID.patch = flow_table[inx].patchID;
			flow_table[pch].stream_inx = inx;
		}
	}
	return;
}

int CreatFlowTable::sort_flow_table(int *patch_pch, double *dem, int *patch, int maxr, int maxc)
{
	// local variable declarations
	int i, j;
	flow_struct temp_entry;

	cout << "\nSTART SORTING FLOW TABLE:\n" << endl;

	// local variable declarations
	int inx{}, iny{};
	int r{}, c{}, pch{};
	double tmp_int{};
	double tmp_dou{};

	num_patches = 0;

	//Before allocating memorys, the patches are initialized with well order
	double *dem_real = new double[maxr*maxc]{};
	int *patch_real = new int[maxr*maxc]{};//patchID
	int *pch_real = new int[maxr*maxc]{};//defines the flow_table list order in 1:num_patches

	cout << "\t initializing::" << endl;
	//initializing
	for (inx = 0; inx < maxr*maxc; inx++) {

		if (patch[inx] > 0) {
			num_patches++;
			dem_real[iny] = dem[inx];
			patch_real[iny] = patch[inx];
			pch_real[iny] = num_patches;//it starts from 1
			iny++;
		}
	}

	cout << "\t sorting:: wait in this process" << endl;
	//sorting
	for (inx = 0; inx < num_patches; inx++) {
		for (iny = num_patches - 1; iny != inx; iny--) {

			if (dem_real[iny] > dem_real[iny - 1]) {

				//exchange dem
				tmp_dou = dem_real[iny];
				dem_real[iny] = dem_real[iny - 1];
				dem_real[iny - 1] = tmp_dou;

				//exchange patch_inx
				tmp_int = patch_real[iny];
				patch_real[iny] = patch_real[iny - 1];
				patch_real[iny - 1] = tmp_int;

			}
		}
		//recording the progress
		if(inx%10000==0)
		cout <<"\t   " <<inx/10000 <<" w"<< endl;
	}
	cout << "\t delivering::it might cost much time with no indicators" << endl;

	//delivering
	for (inx = 0; inx < maxr*maxc; inx++) {

		if (patch[inx] > 0) {

			//seaching and macthing
			for (iny = 0; iny != num_patches; iny++) {
				if (patch[inx] == patch_real[iny]) {
					patch_pch[inx] = iny + 1;//it starts from 1
					break;
				}
			}
		}
	}

	cout << "END OF SORTING FLOW TABLE PATCHES IN AN DESENDING ORDER" << endl;
	return num_patches;
}

void CreatFlowTable::zero_flow_table(flow_struct *flow_table, int maxr, int maxc)
{
	int r, c;
	int inx;

	for (r = 0; r < maxr; r++) {
		for (c = 0; c < maxc; c++) {
			inx = r * maxc + c;
			flow_table[inx].patchID = -99999;
			flow_table[inx].ID_order = -99999;
			flow_table[inx].x = 0.0;
			flow_table[inx].y = 0.0;
			flow_table[inx].z = 0.0;
			flow_table[inx].land = 0;
			flow_table[inx].flna = 0.0;
			flow_table[inx].total_gamma = 0.0;
			flow_table[inx].area = 0;
			flow_table[inx].acc_area = 0.0;
			flow_table[inx].path_length = 0.0;
			flow_table[inx].num_adjacent = 0;
			flow_table[inx].internal_slope = 0;
			flow_table[inx].distance_to_stream = 2000000;
			flow_table[inx].elevation_above_stream = 2000000;
			flow_table[inx].number = 0;
			//flow_table[inx].adj_list = NULL;
			//flow_table[inx].adj_ptr = NULL;
		}
	}
	return;
}

void CreatFlowTable::print_topo_indices(int num_patches, struct flow_struct * flow_table, int *patch,
	CReadCommandline *commandline, double cell, char *input_prefix, double width, int maxr, int maxc, int algorithm_flag)
{
	int pch;
	int inx;
	int pch_now;

	char algorithm[6][10]{ "D8", "MD8", "D_inf", "RMD_inf", "MFD_md", "MD_inf" };

	cout << "\n START PRINTING TOPO INDICES:: SORTED FLOW TABLE TO MATRIX" << endl;

	char pch_dat[256], pch_txt[256];
	char acc_dat[256], acc_txt[256];
	char dts_dat[256], dts_txt[256];
	char sts_dat[256], sts_txt[256];
	char dtr_dat[256], dtr_txt[256];
	char str_dat[256], str_txt[256];
	char tga_dat[256], tga_txt[256];
	char twi_dat[256], twi_txt[256];

	fstream PCH_dat;
	fstream PCH_txt;

	fstream ACC_dat;
	fstream ACC_txt;

	fstream DTS_dat;
	fstream DTS_txt;

	fstream STS_dat;
	fstream STS_txt;

	fstream DTR_dat;
	fstream DTR_txt;

	fstream STR_dat;
	fstream STR_txt;

	fstream TGA_dat;
	fstream TGA_txt;

	fstream TWI_dat;
	fstream TWI_txt;



	if (commandline->pch == true) {

		strcpy(pch_dat, input_prefix);
		strcat(pch_dat, "_pch");
		strcat(pch_dat, ".dat");

		strcpy(pch_txt, input_prefix);
		strcat(pch_txt, "_pch");
		strcat(pch_txt, ".txt");

		PCH_dat.open(pch_dat, ios::out);
		PCH_txt.open(pch_txt, ios::out);

		if (!PCH_dat.is_open() && !PCH_txt.is_open()) {
			cout << "\n opening PCH output file failed" << endl;
			exit(1);
		}
	}



	if (commandline->acc == true) {

		strcpy(acc_dat, input_prefix);
		strcat(acc_dat, "_acc_");
		strcat(acc_dat, algorithm[algorithm_flag]);
		strcat(acc_dat, ".dat");

		strcpy(acc_txt, input_prefix);
		strcat(acc_txt, "_acc_");
		strcat(acc_txt, algorithm[algorithm_flag]);
		strcat(acc_txt, ".txt");

		ACC_dat.open(acc_dat, ios::out);
		ACC_txt.open(acc_txt, ios::out);

		if (!ACC_dat.is_open() && !ACC_txt.is_open()) {
			cout << "\n opening AccArea output file failed" << endl;
			exit(1);
		}

	}


	if (commandline->dts == true) {

		strcpy(dts_dat, input_prefix);
		strcat(dts_dat, "_dts_");
		strcat(dts_dat, algorithm[algorithm_flag]);
		strcat(dts_dat, ".dat");

		strcpy(dts_txt, input_prefix);
		strcat(dts_txt, "_dts_");
		strcat(dts_txt, algorithm[algorithm_flag]);
		strcat(dts_txt, ".txt");

		DTS_dat.open(dts_dat, ios::out);
		DTS_txt.open(dts_txt, ios::out);

		if (!DTS_dat.is_open() && !DTS_txt.is_open()) {
			cout << "\n opening DTS output file failed" << endl;
			exit(1);
		}

	}

	if (commandline->sts == true) {

		strcpy(sts_dat, input_prefix);
		strcat(sts_dat, "_sts_");
		strcat(sts_dat, algorithm[algorithm_flag]);
		strcat(sts_dat, ".dat");

		strcpy(sts_txt, input_prefix);
		strcat(sts_txt, "_sts_");
		strcat(sts_txt, algorithm[algorithm_flag]);
		strcat(sts_txt, ".txt");

		STS_dat.open(sts_dat, ios::out);
		STS_txt.open(sts_txt, ios::out);

		if (!STS_dat.is_open() && !STS_txt.is_open()) {
			cout << "\n opening STS output file failed" << endl;
			exit(1);
		}

	}


	if (commandline->dtr == true) {

		strcpy(dtr_dat, input_prefix);
		strcat(dtr_dat, "_dtr_");
		strcat(dtr_dat, algorithm[algorithm_flag]);
		strcat(dtr_dat, ".dat");

		strcpy(dtr_txt, input_prefix);
		strcat(dtr_txt, "_dtr_");
		strcat(dtr_txt, algorithm[algorithm_flag]);
		strcat(dtr_txt, ".txt");

		DTR_dat.open(dtr_dat, ios::out);
		DTR_txt.open(dtr_txt, ios::out);

		if (!DTR_dat.is_open() && !DTR_txt.is_open()) {
			cout << "\n opening DTR output file failed" << endl;
			exit(1);
		}

	}

	if (commandline->str == true) {

		strcpy(str_dat, input_prefix);
		strcat(str_dat, "_str_");
		strcat(str_dat, algorithm[algorithm_flag]);
		strcat(str_dat, ".dat");

		strcpy(str_txt, input_prefix);
		strcat(str_txt, "_str_");
		strcat(str_txt, algorithm[algorithm_flag]);
		strcat(str_txt, ".txt");

		STR_dat.open(str_dat, ios::out);
		STR_txt.open(str_txt, ios::out);

		if (!STR_dat.is_open() && !STR_txt.is_open()) {
			cout << "\n opening STR output file failed" << endl;
			exit(1);
		}

	}


	if (commandline->tga == true) {

		strcpy(tga_dat, input_prefix);
		strcat(tga_dat, "_tga_");
		strcat(tga_dat, algorithm[algorithm_flag]);
		strcat(tga_dat, ".dat");

		strcpy(tga_txt, input_prefix);
		strcat(tga_txt, "_tga_");
		strcat(tga_txt, algorithm[algorithm_flag]);
		strcat(tga_txt, ".txt");

		TGA_dat.open(tga_dat, ios::out);
		TGA_txt.open(tga_txt, ios::out);

		if (!TGA_dat.is_open() && !TGA_txt.is_open()) {
			cout << "\n opening TGA output file failed" << endl;
			exit(1);
		}
	}

	if (commandline->twi == true) {

		strcpy(twi_dat, input_prefix);
		strcat(twi_dat, "_twi_");
		strcat(twi_dat, algorithm[algorithm_flag]);
		strcat(twi_dat, ".dat");

		strcpy(twi_txt, input_prefix);
		strcat(twi_txt, "_twi_");
		strcat(twi_txt, algorithm[algorithm_flag]);
		strcat(twi_txt, ".txt");

		TWI_dat.open(twi_dat, ios::out);
		TWI_txt.open(twi_txt, ios::out);

		if (!TWI_dat.is_open() && !TWI_txt.is_open()) {
			cout << "\n opening TWI output file failed" << endl;
			exit(1);
		}
	}

	//end of four kinds of topo indicies file creating 


	if (commandline->dts == true || commandline->sts == true || commandline->dtr == true || commandline->str == true) {
		//here we try to sort out the flow steps to stream and distance to stream
		const int max_num = 5000;
		adj_struct *ptr;
		for (int pchA = 1; pchA <= num_patches; pchA++)
		{
			flow_table[pchA].Dts = flow_table[pchA].Sts = 0;
			flow_table[pchA].Dtr = flow_table[pchA].Str = 0;

			//start searching downslope cells
			//processing DIS AND STE at cell level
			ptr = flow_table[pchA].adj_list;
			for (int neigh = 1; neigh <= flow_table[pchA].num_adjacent; neigh++) {

				ptr->Dis = ptr->Ste = 0;
				if (ptr->gamma > 0) {

					ptr->Ste = ptr->gamma;

					for (int pchB = pchA + 1; pchB <= num_patches; pchB++) {

						if (flow_table[pchB].patchID == ptr->patchID) {

							if (pow(flow_table[pchA].x - flow_table[pchB].x, 2.0) < 0.01 || pow(flow_table[pchA].y - flow_table[pchB].y, 2.0) < 0.01)
								ptr->Dis = ptr->gamma*cell;
							else
								ptr->Dis = ptr->gamma*cell*pow(2, 1 / 2.0);
							break;
						}
					}
				}
				ptr = ptr->next;
			}
		}



		//processing DTS AND STS
		if (commandline->dts == true || commandline->sts == true) {
			for (int pchA = num_patches; pchA != 0; pchA--)
			{

				//start searching downslope cells
				if (flow_table[pchA].land != 1) {
					ptr = flow_table[pchA].adj_list;
					for (int neigh = 1; neigh <= flow_table[pchA].num_adjacent; neigh++) {

						if (ptr->gamma > 0) {

							flow_table[pchA].Dts += ptr->Dis;
							flow_table[pchA].Sts += ptr->Ste;
							for (int pchB = pchA + 1; pchB <= num_patches; pchB++) {

								//matching
								if (flow_table[pchB].patchID == ptr->patchID) {

									flow_table[pchA].Dts += flow_table[pchB].Dts* ptr->gamma;
									flow_table[pchA].Sts += flow_table[pchB].Sts* ptr->gamma;
									break;
								}

							}

						}
						ptr = ptr->next;
					}
				}
			}
		}


		//processing DTR AND STR
		if (commandline->dtr == true || commandline->str == true) {
			for (int pchA = 1; pchA <= num_patches; pchA++)
			{

				//start searching downslope cells
				int ok = 0;
				//if (flow_table[pchA].land != 1) {
				if (ok != 1) {
					ptr = flow_table[pchA].adj_list;
					for (int neigh = 1; neigh <= flow_table[pchA].num_adjacent; neigh++) {

						if (ptr->gamma > 0) {

							for (int pchB = pchA + 1; pchB <= num_patches; pchB++) {

								//matching
								if (flow_table[pchB].patchID == ptr->patchID) {

									flow_table[pchB].Dtr += flow_table[pchA].Dtr*ptr->gamma;
									flow_table[pchB].Dtr += ptr->Dis;

									flow_table[pchB].Str += flow_table[pchA].Str* ptr->gamma;
									flow_table[pchB].Str += ptr->Ste;

									break;

								}

							}

						}
						ptr = ptr->next;
					}
				}
			}
		}

	}//END OF DTS && STS


	//processing DTR AND STR
	if (commandline->twi == true) {
		for (int pchA = 1; pchA <= num_patches; pchA++)
		{
			//start searching downslope cells
			flow_table[pchA].TWI = log(flow_table[pchA].acc_area*cell / flow_table[pchA].total_gamma);
		}
	}



	//THIS STEP ONLY DILIVER VALUES
	int *patch_rec = new int[maxc*maxr];
	for (int r = 0; r < maxr; r++)
		for (int c = 0; c < maxc; c++)
		{
			inx = r * maxc + c;
			patch_rec[inx] = patch[inx];//we changed here
		}

	for (int r = 0; r < maxr; r++) {
		for (int c = 0; c < maxc; c++) {
			//cout << c << "\t" << r << endl;
			inx = r * maxc + c;//rectified
			if (patch[inx] == -9999) {

				if (commandline->acc == true)
					ACC_dat << setw(14) << "NA";

				if (commandline->dts == true)
					DTS_dat << setw(14) << "NA";

				if (commandline->sts == true)
					STS_dat << setw(14) << "NA";

				if (commandline->dtr == true)
					DTR_dat << setw(14) << "NA";

				if (commandline->str == true)
					STR_dat << setw(14) << "NA";

				if (commandline->tga == true)
					TGA_dat << setw(14) << "NA";

				if (commandline->twi == true)
					TWI_dat << setw(14) << "NA";

				if(commandline->pch==true)
					PCH_dat << setw(14) << "NA";

			}
			else
			{
				for (pch = 1; pch <= num_patches; pch++)//SEARCH VALUE
				{
					if (flow_table[pch].patchID == patch[inx])
					{
						pch_now = pch;

						if (commandline->pch == true){
							PCH_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << pch-1;//start from 0
							PCH_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << pch-1 << endl;//start from 0
						}


						if (commandline->acc == true) {
							ACC_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].acc_area;
							ACC_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].acc_area << endl;

						}

						if (commandline->dts == true) {
							DTS_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Dts;
							DTS_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Dts << endl;

						}

						if (commandline->sts == true) {
							STS_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Sts;
							STS_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Sts << endl;
						}

						if (commandline->dtr == true) {
							DTR_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Dtr;
							DTR_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Dtr << endl;
						}

						if (commandline->str == true) {
							STR_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Str;
							STR_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].Str << endl;
						}

						if (commandline->tga == true) {
							TGA_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].total_gamma;
							TGA_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].total_gamma << endl;

						}

						if (commandline->twi == true) {

							TWI_dat << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].TWI;
							TWI_txt << setw(14) << setiosflags(ios::fixed) << setprecision(2) << flow_table[pch].TWI << endl;
						}


						break;
					}
				}
			}
		}//END OF A ROW
		if (commandline->acc == true) ACC_dat << endl;
		if (commandline->dts == true) DTS_dat << endl;
		if (commandline->sts == true) STS_dat << endl;
		if (commandline->dtr == true) DTR_dat << endl;
		if (commandline->str == true) STR_dat << endl;
		if (commandline->tga == true) TGA_dat << endl;
		if (commandline->twi == true) TWI_dat << endl;
		if (commandline->pch == true) PCH_dat << endl;

	}//END OF ALL 

	//Closing files
	delete patch_rec;

	if (commandline->pch == true) {
		PCH_dat.close();
		PCH_txt.close();
	}

	if (commandline->acc == true) {
		ACC_dat.close();
		ACC_txt.close();
	}
	if (commandline->dts == true) {
		DTS_dat.close();
		DTS_txt.close();
	}
	if (commandline->sts == true) {
		STS_dat.close();
		STS_txt.close();
	}
	if (commandline->dtr == true) {
		DTR_dat.close();
		DTR_txt.close();
	}
	if (commandline->str == true) {
		STR_dat.close();
		STR_txt.close();
	}
	if (commandline->tga == true) {
		TGA_dat.close();
		TGA_txt.close();
	}

	if (commandline->twi == true) {
		TWI_dat.close();
		TWI_txt.close();
	}

	cout << "\n END OF PRINTING TOPO INDICES" << endl;
	return;
}

void	CreatFlowTable::print_flow_table(int num_patches, flow_struct *flow_table, int sc_flag, double cell, char *input_prefix, double width, int algorithm_flag)
{
	int i, j, cnt;
	adj_struct *adj_ptr;

	double mult, tmp, unidirectional;
	static int count = 0;
	char name[256];
	char algorithm[6][10]{ "D8", "MD8", "D_inf", "RMD_inf", "MFD_md","MD_inf" };
	ofstream FlowTable;

	strcpy(name, input_prefix);
	strcat(name, "_flow_table_");
	strcat(name, algorithm[algorithm_flag]);
	strcat(name, ".dat");


	FlowTable.open(name, ios::out);

	if (!FlowTable) {
		cout << "\n opening flow_table output file failed" << endl;
		exit(1);
	}

	FlowTable << num_patches << endl;
	for (i = 1; i <= num_patches; i++) {
		count = 0;

		// this is a temporary patch, so that streams immediately produces outflow
		// below are for stream patch
		if ((flow_table[i].land == 1) && (sc_flag > 0)) {
			if (flow_table[i].land == 1) {
				if (flow_table[i].area > 1) {
					std::cout << "flow_table[i].area is " << flow_table[i].area << std::endl;
					getchar();
				}

				//mult = ( double )( flow_table[i].K * flow_table[i].m_par * sqrt(flow_table[i].area) ); 	
				//mult = ( double )(0.12*2* sqrt(flow_table[i].area) );
				mult = 1.;
				if (sc_flag == 2)
					tmp = (double)(rand() / (pow(2.0, 15.0) - 1));
				else
					tmp = flow_table[i].internal_slope;

				if (sc_flag == 3) {
					adj_ptr = flow_table[i].adj_list;
					flow_table[i].internal_slope = 0.0;
					cnt = 0;
					for (j = 1; j <= flow_table[i].num_adjacent; j++) {
						if (adj_ptr->gamma <= 0)
							flow_table[i].internal_slope += adj_ptr->slope;
						cnt += 1;
						adj_ptr = adj_ptr->next;
					}
					flow_table[i].internal_slope = flow_table[i].internal_slope / cnt;
					tmp = (double)(-1.0*flow_table[i].internal_slope);
				}
				flow_table[i].total_gamma = (double)(mult * tmp);
			}
		}

		unidirectional = 0.001;
		adj_ptr = flow_table[i].adj_list;
		for (j = 1; j <= flow_table[i].num_adjacent; j++) {
			if (adj_ptr->gamma == 0.)
				count++;

			//Below if for judge unidirectional flow
			if (adj_ptr->gamma > unidirectional)
			{
				unidirectional = adj_ptr->gamma;
			}
			adj_ptr = adj_ptr->next;
		}

		if (flow_table[i].total_gamma < 0.001) flow_table[i].total_gamma = 0.001;


		FlowTable << fixed << right << setw(8) << flow_table[i].patchID
			<< setw(12) << setprecision(1) << flow_table[i].x
			<< setw(12) << setprecision(1) << flow_table[i].y
			<< setw(12) << setprecision(3) << flow_table[i].z
			<< setw(10) << setprecision(2) << flow_table[i].acc_area
			<< setw(10) << flow_table[i].land
			<< setw(14) << setprecision(10) << flow_table[i].total_gamma
			<< setw(4) << flow_table[i].num_adjacent - count << endl;
		
		if (flow_table[i].total_gamma < 0.001) cout << flow_table[i].patchID<<"\t"<< flow_table[i].total_gamma << endl;

		/*
		fprintf_s(outfile,"\n %-8d %-8.1f %-8.1f %-7.1f %-10.2f %-4d %-12.5f %-4d",
		flow_table[i].patchID,
		flow_table[i].x,
		flow_table[i].y,
		flow_table[i].z,
		flow_table[i].acc_area,
		flow_table[i].land,
		flow_table[i].total_gamma,
		1);
		*/

		adj_ptr = flow_table[i].adj_list;
		for (j = 1; j <= flow_table[i].num_adjacent; j++) {

			if (adj_ptr->gamma > 0.)//inx-1 equals to patchorder starts from 0
				FlowTable << setw(13) << adj_ptr->inx-1<< setw(13) << adj_ptr->patchID << setw(13) << setprecision(8) << adj_ptr->gamma << endl;
			adj_ptr = adj_ptr->next;
		}

		if (flow_table[i].land == 2) {
			FlowTable << setw(16) << flow_table[i].stream_ID.patch << setw(1) << width << endl;
		}
	}

	FlowTable.close();

	return;
}

void	CreatFlowTable::print_drain_stats(int num_patches, struct flow_struct *flow_table, char *input_prefix)
{
	int i;
	ofstream outfile;
	char name[256];

	std::cout << "\n Printing Drainage Statistics " << std::endl;

	strcpy(name, input_prefix);
	strcat(name, "_stats.dat");
	outfile.open(name, ios::out);
	if (!outfile) {
		std::cout << "Error opening flow_table output file\n" << std::endl;
		exit(1);
	}

	outfile << setw(8) << num_patches << endl;
	for (i = 1; i <= num_patches; i++) {

		outfile << setw(6) << flow_table[i].patchID
			<< flow_table[i].x
			<< flow_table[i].y
			<< flow_table[i].z
			<< flow_table[i].patchID
			<< flow_table[i].acc_area
			<< flow_table[i].land
			<< flow_table[i].total_gamma
			<< flow_table[i].slope
			<< flow_table[i].area
			<< flow_table[i].num_adjacent << endl;
	}

	//fprintf_s(outfile,"\n %6d %6d %6d %6.1f %6.1f %6.1f %10f %4d %15.8f %10.6f %9d %4d",  

	outfile.close();

	return;
}

#endif