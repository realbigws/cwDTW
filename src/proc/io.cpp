#include "io.h"

bool g::io::ReadATCG(const char* name, std::vector<char>& genomes)
{
    std::ifstream in(name);
    if(!in.good()) {
        return false;
    }
    
    //-> skip first header
    std::string buf;
    if(!getline(in, buf)){
        return false;
    }
    
    while(in.good()){
        char item;
        in>>item;
        if(in.fail()){
            break;
        }
        genomes.push_back(item);
    }
    in.close();
	
    return true;
}

bool g::io::WriteATCG(const char* name, const std::vector<char>& genomes)
{
	std::ofstream out(name);
    if(!out.good()) {
        return false;
    }
    
    for(size_t i = 0; i < genomes.size(); i++){
		out<<genomes[i]<<std::endl;
    }
    
    out.close();
    
    return true;
}

bool g::io::ReadSignalSequence(const char* name, std::vector<double>& signals)
{
    std::ifstream in(name);
    if(!in.good()) {
        return false;
    }

    while(in.good()){
        double item;
        in>>item;
        if(in.fail()){
            break;
        }
        signals.push_back(item);
    }
    in.close();
	
    return true;
}

bool g::io::ReadSignalSequence_int(const char* name, std::vector<int>& signals)
{
    std::ifstream in(name);
    if(!in.good()) {
        return false;
    }

    while(in.good()){
        int item;
        in>>item;
        if(in.fail()){
            break;
        }
        signals.push_back(item);
    }
    in.close();

    return true;
}

bool g::io::WriteSignalSequence(const char* name, const std::vector<int>& signals)
{
	std::ofstream out(name);
    if(!out.good()) {
        return false;
    }
    
    for(size_t i = 0; i < signals.size(); i++){
		out<<signals[i]<<std::endl;
    }
    
    out.close();
    
    return true;
}
