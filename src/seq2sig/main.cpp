#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <proc/io.h>
#include "util/exception.h"
#include "5mer/5mer_index.h"

bool Genomes2SignalSequence(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<int>& signals, int scale, int FIVE_or_SIX)
{
	size_t bound;
	if(FIVE_or_SIX==0) //-> 5mer model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-5;//genomes.size()%5;
		for(size_t i = 0; i < bound; i++){
			double sigval = g::Mer2Signal::AvgSignalAt_5mer(index[i]);
			sigval = 5.7*sigval+14;
			for(int c = scale; c--;){
				signals.push_back(sigval);
			}
		}
	}
	else
	{
		g::Mer2Signal::Genome2Index_6mer(genomes, index);
		bound = genomes.size()-6;//genomes.size()%5;
		for(size_t i = 0; i < bound; i++){
			double sigval = g::Mer2Signal::AvgSignalAt_6mer(index[i]);
			sigval = 5.7*sigval+14;
			for(int c = scale; c--;){
				signals.push_back(sigval);
			}
		}
	}

	//---- tail five_mer ------//
//	for(size_t i = bound; i < genomes.size(); i++){
//		for(int c = scale; c--;){
//			signals.push_back(100);
//		}
//	}
}

//---------------------- main ---------------------//
int main(int argc, char **argv)
{
	std::string input="";
	std::string output="";

	struct options opts;
	opts.scale=1;
	opts.kmer=0;
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	input=opts.input;
	output=opts.output;
	if(input=="" || output=="")
	{
		fprintf(stderr,"input or output is NULL \n");
		return -1;
	}

	std::vector<char> genomes;
	std::vector<int> index;
	std::vector<int> signals;
	g::io::ReadATCG(opts.input, genomes);
	Genomes2SignalSequence(genomes, index, signals, opts.scale, opts.kmer);
	g::io::WriteSignalSequence_int(opts.output, signals);

	return 0;
}

