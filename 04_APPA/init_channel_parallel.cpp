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

void init_channel_parallel(struct patch_object* patch, struct status_information* INF){


	//prepare channel information
	for (int pch = 0; pch != INF->patch_num; pch++) {

			patch[pch].channel_state = 0;//as undefined
			patch[pch].channel_inx = 0;
			patch[pch].layer_thread = 0;
			patch[pch].re_channel_acc = patch[pch].channel_acc;//init as original channel_acc
			patch[pch].lock_state = 0;

	}

	return;
}