#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/cwdtw.h"
#include "util/exception.h"
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



//----------- write alignment to file -------------//__2017.10.15__//(Sheng modified)
void WriteSequenceAlignment(const char* output, 
	const std::vector<double>& reference_orig, const std::vector<double>& peer_orig,
	const std::vector<double>& reference, const std::vector<double>& peer,
	vector<pair<long,long> >& alignment, int swap)
{
	vector <std::string> tmp_rec;
	double diff;	
	for(long i = 0; i < alignment.size(); i++)
	{
		//----- output to string ----//
		std::ostringstream o;
		diff = std::fabs(reference[alignment[i].first]-peer[alignment[i].second]);
		if(swap==1)
		{
			o<<setw(10)<<alignment[i].second+1<<" "<<setw(10)<<alignment[i].first+1<<" | ";
			o<<setw(15)<<reference_orig[alignment[i].second]<<" "<<setw(15)<<peer_orig[alignment[i].first]<<" | ";
			o<<setw(15)<<peer[alignment[i].second]<<" "<<setw(15)<<reference[alignment[i].first];
		}
		else
		{
			o<<setw(10)<<alignment[i].first+1<<" "<<setw(10)<<alignment[i].second+1<<" | ";
			o<<setw(15)<<reference_orig[alignment[i].first]<<" "<<setw(15)<<peer_orig[alignment[i].second]<<" | ";
			o<<setw(15)<<reference[alignment[i].first]<<" "<<setw(15)<<peer[alignment[i].second];
		}
		o<<"          diff:"<<setw(15)<<diff;
		//----- record string -----//
		std::string s=o.str();
		tmp_rec.push_back(s);
	}
	//----- output to file ------//
	FILE *fp=fopen(output,"wb");
	for(long i=0;i<(long)tmp_rec.size();i++)fprintf(fp,"%s\n",tmp_rec[i].c_str());
	fclose(fp);
}

//----------------- main -------------------//
int main(int argc, char **argv)
{
	struct options opts;
	opts.radius  = 50;
	opts.level   = 3;
	opts.scale0  = sqrt(2);
	opts.verbose = 0;       //-> [0] no verbose; 1 verbose
	opts.test    = 0;       //-> [0] not use test mode; 1 equal_ave, 2 peak_ave, 3 Fast_DTW
	opts.mode    = 0;       //-> [0] block bound; 1 diagonol bound
	

	//----- parse arguments -----//
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	std::string input1=opts.input;
	std::string input2=opts.peer;
	std::string output=opts.output;
	if(input1=="" || input2=="")
	{
		fprintf(stderr,"input1 or input2 is NULL \n");
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
	std::vector<double> reference_orig=reference;
	std::string genom_name_orig=opts.input;
	std::string genom_name;
	getBaseName(genom_name_orig,genom_name,'/','.');

	//========================================//
	//----- 2. read nanopore raw signal -----//
	std::vector<double> peer;
	if(!g::io::ReadSignalSequence(opts.peer,peer)){
		EX_TRACE("Cannot open %s.\n", opts.peer);
		return -1;
	}
	std::vector<double> peer_orig=peer;
	std::string signal_name_orig=opts.peer;
	std::string signal_name;
	getBaseName(signal_name_orig,signal_name,'/','.');


	//----- length check ------//
	int swap=0;
	if(reference.size()>peer.size())
	{
		std::vector<double> tmp=peer;
		peer=reference;
		reference=tmp;
		swap=1;
	}

	//==================================================//
	//------3. process initial input signals ----------//

	//----- 3.1 Zscore normaliza on both signals -----//
	g::proc::ZScoreNormalize(reference);
	g::proc::ZScoreNormalize(peer);

	//----- 3.2 calculate length ratio between input signals -----//	
	double alpha = (double)peer.size()/reference.size();


	//--------------- Just For Test (Fast_DTW) -------------------//
	if(opts.test==3)
	{
		//------ 5.1 generate initial alignment via FastDTW ------//
		vector<pair<long,long> > cosali;
		double tdiff;
		long max_level=100000; //here we fix maximal level 
		g::cwdtw::FastDTW(reference, peer, cosali, opts.radius, max_level,opts.mode, &tdiff);
		vector<pair<long,long> > bound;
		g::cwdtw::BoundGeneration(cosali, opts.radius, bound, opts.mode);
		//------ 5.2 generate final alignment via cDTW ------//
		std::vector<std::pair<long,long> > alignment;
		tdiff = g::proc::BoundDynamicTimeWarping(reference, peer, bound, alignment);
		fprintf(stderr,"%s %s %lf %d %lf\n",signal_name.c_str(),genom_name.c_str(),tdiff,alignment.size(),tdiff/alignment.size());

		//=================================================//
		//------ 6. output final alignment to file -------//
		if(output!="")
			WriteSequenceAlignment(opts.output, reference_orig, peer_orig, reference, peer, alignment, swap);
			
		//----- exit -----//	
		return 0;
	}


	//====================================================//
	//----- 4. continous wavelet transform --------------//
	std::vector<std::vector<double> > rcwt, pcwt;

	if(opts.verbose ==1){
		EX_TRACE("CWT Analysis...\n");
	}
	
	long npyr = opts.level;          // default: 3
	double scale0 = opts.scale0;	// default: sqrt(2)
	double dscale = 1;              // default: 1

	g::cwdtw::CWTAnalysis(reference, rcwt, scale0, dscale, npyr);	
	g::cwdtw::CWTAnalysis(peer, pcwt, scale0*alpha, dscale, npyr);

	//------ 4.1 Zscore normaliza on both CWT signals -----//	
	//if multiscale is used, pyr logical should be added.
	for(long i = 0; i < npyr; i++){
		g::proc::ZScoreNormalize(rcwt[i]);
		g::proc::ZScoreNormalize(pcwt[i]);
	}


	//============================================//
	//------ 5. multi-level WaveletDTW ----------//
	std::vector<std::pair<long,long> > cosali;
	double tdiff;

	if(opts.verbose ==1){
		EX_TRACE("Coarse Alignment...\n");
	}
	g::cwdtw::MultiLevel_WaveletDTW(reference, peer, rcwt, pcwt, cosali, opts.radius, opts.test, opts.mode, &tdiff);
	if(opts.verbose ==1){
		EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff/cosali.size());
	}

	//------ 5.1 generate final boundary -------//
	std::vector<std::pair<long,long> > bound;
	g::cwdtw::BoundGeneration(cosali, opts.radius, bound, opts.mode);

	//------ 5.2 generate final alignment via cDTW ------//
	std::vector<std::pair<long,long> > alignment;
	tdiff = g::proc::BoundDynamicTimeWarping(reference, peer, bound, alignment);
	fprintf(stderr,"%s %s %lf %d %lf\n",signal_name.c_str(),genom_name.c_str(),tdiff,alignment.size(),tdiff/alignment.size());


	//=================================================//
	//------ 6. output final alignment to file -------//
	if(output!="")
		WriteSequenceAlignment(opts.output, reference_orig, peer_orig, reference, peer, alignment, swap);

	//----- exit -----//	
	return 0;
}

