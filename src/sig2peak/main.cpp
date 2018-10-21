#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/cwdtw.h"
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



//----------------- main -------------------//
int main(int argc, char **argv)
{
	struct options opts;
	opts.scale0 = sqrt(2);
	opts.ZorNOT = 0;    //0 for NOT perform Z-normalize
	opts.posout = 0;    //0 for NOT output the position of each peak

	//----- parse arguments -----//
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	std::string input="";
	std::string output="";
	input=opts.input;
	output=opts.output;
	if(input=="" || output=="")
	{
		fprintf(stderr,"input or output is NULL \n");
		return -1;
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


//printf("read done \n");

	//==================================================//
	//----- 2. process initial input signals ----------//
	if(opts.ZorNOT==1){
		g::proc::ZScoreNormalize(reference);
	}

	//====================================================//
	//----- 3. continous wavelet transform --------------//
	std::vector<std::vector<double> > rcwt;
	long npyr = 1;                   // default: 1
	double scale0 = opts.scale0;	// default: sqrt(2)
	double dscale = 1;              // default: 1
	g::cwdtw::CWTAnalysis(reference, rcwt, scale0, dscale, npyr);	

//printf("proc CWT done \n");

	//----- 4. Zscore normaliza on both CWT signals -----//	
	//if multiscale is used, pyr logical should be added.
	if(opts.ZorNOT==1){
		g::proc::ZScoreNormalize(rcwt[0]);
	}
	std::vector<std::pair<long, double> > sigpeaks;
	g::proc::PeakPick(rcwt[0], sigpeaks);

//printf("proc PeakPick done \n");

	//=================================================//
	//----- 5. output final alignment to file -------//
	FILE *fp=fopen(output.c_str(),"wb");
	for(long i=0;i<sigpeaks.size();i++){
		if(opts.posout==0){
			fprintf(fp,"%lf\n",sigpeaks[i].second);
		}
		else{
			fprintf(fp,"%lf %d\n",sigpeaks[i].second,sigpeaks[i].first);
		}
	}
	fclose(fp);

//printf("write done \n");

	//----- exit -----//	
	return 0;
}

