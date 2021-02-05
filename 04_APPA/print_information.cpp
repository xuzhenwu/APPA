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

void print_information(struct status_information *INF, int available_threads)
{
	// compute ESR and TMSR
	double pch = 1. * INF->channel_num / INF->patch_num;
	double phs = 1. * (INF->patch_num - INF->channel_num) / INF->patch_num;
	INF->ESR  = 1 / (phs / INF->ESRhs + pch / INF->ESRch);
	INF->TMSR = 1 / (phs / INF->TMSRhs + pch / min(INF->TMSRch, 1. * available_threads)); // modify for TMSRch > available_threads 

	cout << "Threads:" << available_threads 
		<< "\tESR:" << fixed << setprecision(2) << INF->ESR 
		<< "\tESRhs:"<< fixed << setprecision(2) << INF->ESRhs 
		<< "\tESRch:" << fixed << setprecision(2) << INF->ESRch
		<< "\tTMSR:" << fixed << setprecision(2) << INF->TMSR
		<< "\tTMSRhs:" << fixed << setprecision(2) << INF->TMSRhs
		<< "\tTMSRch:" << fixed << setprecision(2) << min(INF->TMSRch, 1. * available_threads)
		<< "\n" << endl;
}

