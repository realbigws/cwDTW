#ifndef FIVEMER_INDEX_H__
#define FIVEMER_INDEX_H__

#include <vector>
#include <string>

namespace g{

class Mer2Signal
{
private:
	//---- in-house made 5-mer pore model ----//
	/* avg, var, len, #, #, #, # (defined in 5mer_index_table.rc) */
	//const static double index_table_self[1024][6];   //-> (obsoleted)

	//---- official 5/6-mer pore model ------//
	//-> download from 'https://github.com/nanoporetech/kmer_models/tree/master/r9.4_180mv_450bps_6mer'
	const static double index_table_5mer[1024];
	const static double index_table_6mer[4096];
	//-> download from:
	// 'https://github.com/nanoporetech/kmer_models/blob/master/r9.4_200mv_70bps_5mer_RNA/template_median68pA.model'
	// 'https://github.com/nanoporetech/kmer_models/blob/master/r9.4_180mv_70bps_5mer_RNA/template_median69pA.model'
	const static double rna_table_5mer[1024][2];
	
public:
	//----- Nucleotide_to_Int ---//
	static int Nucleotide_to_Int(char c);

	//----- generate index -------//
	static void Genome2Index_5mer(const std::vector <char> &input, std::vector <int> &index);
	static void Genome2Index_6mer(const std::vector <char> &input, std::vector <int> &index);

	//----- extract content -----//
	//static double AvgSignalAt_SelfModel(int index);   //-> (obsoleted)
	static double AvgSignalAt_5mer(int index);
	static double AvgSignalAt_6mer(int index);
	//-> RNA
	static double  RnaSignalAt_5mer_200mv(int index);
	static double  RnaSignalAt_5mer_180mv(int index);
};

}

#endif

