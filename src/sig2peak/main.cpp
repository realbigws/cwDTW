#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/proc.h"
#include "util/exception.h"
#include "wavelib.h"
#include <malloc.h>
#include <cmath>
#include <iomanip>

using namespace std;
using namespace g::proc;

//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}
//---------------------- utility ------//over



//--------------- continuous wavelet transform (CWT) analysis -----------------//
/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void CWTAnalysis(const std::vector<double>& raw, std::vector<std::vector<double> >& output, double scale0, double dscale, int npyr)
{
	const double* sigs = &raw[0];		//sst_nino3.dat
	cwt_object wt;

	size_t N = raw.size();
	double dt = 1;//2;		//sample rate	>  maybe we should use 2?
// 	npyr =  1; 			// Total Number of scales
	wt = cwt_init("dog", 2.0, N, dt, npyr);	//"morlet", "dog", "paul"
	setCWTScales(wt, scale0, dscale, "pow", 2.0);
	cwt(wt, sigs);

	output.resize(npyr);
	for(size_t k = npyr; k--;){
		int idx = npyr-k-1;
		
		output[idx].resize(raw.size());
		size_t offset = k*raw.size();
		for(size_t i = 0; i < output[idx].size(); i++){
			output[idx][i] = wt->output[i+offset].re;
		}
	}
	cwt_free(wt);
}

//----------------- main -------------------//
int main(int argc, char **argv)
{
	struct options opts;
	opts.scale0 = sqrt(2);
	opts.ZorNOT = 0;    //0 for NOT perform Z-normalize

	//----- parse arguments -----//
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return 0;
	}

	std::string input="";
	std::string output="";
	input=opts.input;
	output=opts.output;
	if(input=="" || output=="")
	{
		fprintf(stderr,"input or output is NULL \n");
		exit(-1);
	}

	//======================= START Procedure ===================================//

	//=========================================//
	//----- 1. read genome translated signal -----//
	std::vector<double> reference;
	if(!g::io::ReadSignalSequence(opts.input,reference)){
		EX_TRACE("Cannot open %s.\n", opts.input);
		return -1;
	}
	std::string genom_name_orig=opts.input;
	std::string genom_name;
	getBaseName(genom_name_orig,genom_name,'/','.');

	//==================================================//
	//----- 2. process initial input signals ----------//
	if(opts.ZorNOT==1){
		g::proc::ZScoreNormalize(reference);
	}

	//====================================================//
	//----- 3. continous wavelet transform --------------//
	std::vector<std::vector<double> > rcwt;
	int npyr = 1;                   // default: 1
	double scale0 = opts.scale0;	// default: sqrt(2)
	double dscale = 1;              // default: 1
	CWTAnalysis(reference, rcwt, scale0, dscale, npyr);	

	//----- 4. Zscore normaliza on both CWT signals -----//	
	//if multiscale is used, pyr logical should be added.
	if(opts.ZorNOT==1){
		g::proc::ZScoreNormalize(rcwt[0]);
	}
	std::vector<std::pair<int, double> > sigpeaks;
	g::proc::PeakPick(rcwt[0], sigpeaks);

	//=================================================//
	//----- 5. output final alignment to file -------//
	FILE *fp=fopen(output.c_str(),"wb");
	for(int i=0;i<sigpeaks.size();i++){
		fprintf(fp,"%lf\n",sigpeaks[i].second);
	}
	fclose(fp);

	//----- exit -----//	
	return 0;
}

