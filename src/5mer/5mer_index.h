#ifndef FIVEMER_INDEX_H__
#define FIVEMER_INDEX_H__

#include <vector>
#include <string>

namespace g{

class Mer2Signal
{
private:
	/* avg, var, len, #, #, #, # (defined in 5mer_index_table.rc) */
	const static double index_table[1024][6];
	
public:
	static int Nucleotide_to_Int(char c);
	static void Genome2Index(const std::vector<char> &input, std::vector <int> &index);

	//----- fivemer to index --------//
	static int FiveMer2Index(const std::vector<char>& fivemer);
	static int FiveMer2Index(char g0, char g1, char g2, char g3, char g4);
	static int FiveMer2Index(const std::string& fivemer);
	static int FiveMer2Index(char const* fivemer);
	static double AvgSignalAt(int index);
};

}

#endif

