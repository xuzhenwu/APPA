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

void init_subbasin_parallel(struct patch_struct *patch,int patch_num) {

	for (int pch = 0; pch != patch_num; pch++) {
		//init re_acc_area
		patch[pch].re_acc_area = patch[pch].acc_area;
		patch[pch].basin_inx = 0;
		patch[pch].sthread= 0;
	}
	return;
}