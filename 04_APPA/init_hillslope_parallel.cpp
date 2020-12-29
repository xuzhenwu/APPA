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

void init_hillslope_parallel(struct patch_struct *patch,struct inf* INF) {

	for (int pch = 0; pch != INF->patch_num; pch++) {
		//init re_acc_area
		patch[pch].re_hill_acc = patch[pch].hill_acc;
		patch[pch].basin_inx = 0;
		patch[pch].sthread= 0;
	}
	return;
}