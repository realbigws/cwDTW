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

bool Genomes2SignalSequence(const std::vector<char>& genomes, std::vector<double>& signals, int scale)
{
	size_t bound = genomes.size()-5;//genomes.size()%5;
	for(size_t i = 0; i < bound; i++){
		int idx = g::Mer2Signal::FiveMer2Index(genomes[i], genomes[i+1], genomes[i+2], genomes[i+3], genomes[i+4]);
		double sigval = g::Mer2Signal::AvgSignalAt(idx);
		
		for(int c = scale; c--;){
			signals.push_back(sigval);
		}
	}
	
	for(size_t i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			signals.push_back(100);
		}
	}
}

int main(int argc, char **argv)
{
	std::string input="";
	std::string output="";

	struct options opts;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}

	input=opts.input;
	output=opts.output;
	if(input=="" || output=="")
	{
		fprintf(stderr,"input or output is NULL \n");
		exit(-1);
	}

	
	EX_TRACE("Transform genomes to signal sequence...\n");
	
	std::vector<char> genomes;
	std::vector<double> signals;
	
	g::io::ReadATCG(opts.input, genomes);
	EX_TRACE("%ld genomes are readed.\n", genomes.size());
	
	Genomes2SignalSequence(genomes, signals, 1);
	
	g::io::WriteSignalSequence(opts.output, signals);
	
// 	genome::Mer2Signal::FiveMer2Index("ATAAA");
	
	return 0;
}

