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

void update_channel_parallel_within_layer(struct patch_struct *patch, int patch_num, int layer, int init_flag) {

	if (init_flag == 0)//renew transient state
		for (int pch = 0; pch != patch_num; pch++) {
			//init re_acc_area

			if (patch[pch].channel_state == 1 )
			{
				patch[pch].channel_state = 0;//as undefined
				patch[pch].channel_inx = 0;
				patch[pch].layer_thread = 0;
				patch[pch].re_channel_acc = patch[pch].re_channel_acc_backup;//retained as backup value
			}
			if (patch[pch].lock_state != 0) {
				patch[pch].re_channel_acc = patch[pch].re_channel_acc_backup; //retained as backup value
			}
			patch[pch].lock_state = 0;//unlocked all pch

		}
	else//as fianl thread
		for (int pch = 0; pch != patch_num; pch++) {
			//init re_acc_area

			if (patch[pch].channel_state == 1 && patch[pch].landID == 1)
			{
				patch[pch].channel_state = 2;//as static
				patch[pch].layer = layer+1;//start from 1
				patch[pch].cthread = (layer+1) * 1000 + patch[pch].layer_thread;//into fianl numbers for cthread
			}
			patch[pch].lock_state = 0;
		}

	
	for (int pch = 0; pch != patch_num; pch++) {

		if (patch[pch].channel_state == 0 && patch[pch].landID == 1)	
		{
			patch[pch].re_channel_acc = 1;
		}
	}
	for (int pch = 0; pch != patch_num; pch++) {

		if (patch[pch].channel_state == 0 && patch[pch].landID == 1)
		{
			for (int neigh = 0; neigh < patch[pch].neigh_num; neigh++) {
			
				patch[patch[pch].neigh_pch[neigh]].re_channel_acc += patch[pch].re_channel_acc*patch[pch].neigh_gamma[neigh];
			
			}
			
		}
	}
	


	return;
}