#ifndef CWDTW_H__
#define CWDTW_H__

#include "io.h"
#include "proc.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "wavelib.h"

using namespace std;

namespace g{
namespace cwdtw{

void CWTAnalysis(
	const std::vector<double>& raw, 
	std::vector<std::vector<double> >& output, 
	double scale0, double dscale, long npyr);

void BoundGeneration(
	std::vector<std::pair<long,long> >& cosali, 
	long neib, std::vector<std::pair<long,long> >& bound, int mode,
	int RENMIN_or_SHENG=0);

void FastDTW(
	std::vector<double>& in1, std::vector<double>& in2,
	std::vector<std::pair<long,long> >& alignment, 
	long radius, long max_level, int mode,
	double* totaldiff = 0);

void MultiLevel_WaveletDTW(
	std::vector<double>& in1, std::vector<double>& in2,
	std::vector<std::vector<double> >& sig1, 
	std::vector<std::vector<double> >& sig2, 
	std::vector<std::pair<long,long> >& alignment, 
	long radius, int test, int mode, 
	double* totaldiff = 0);

}

}


#endif

