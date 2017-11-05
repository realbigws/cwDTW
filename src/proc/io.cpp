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
    std::ostringstream o;
    for(size_t i = 0; i < genomes.size(); i++){
        o<<genomes[i]<<std::endl;
    }
    std::string s=o.str();

    //--- output ----//
    FILE *fp=fopen(name,"wb");
    fprintf(fp,"%s",s.c_str());
    fclose(fp);

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

bool g::io::WriteSignalSequence(const char* name, const std::vector<double>& signals)
{
    std::ostringstream o;
    for(size_t i = 0; i < signals.size(); i++){
        o<<signals[i]<<std::endl;
    }
    std::string s=o.str();

    //--- output ----//
    FILE *fp=fopen(name,"wb");
    fprintf(fp,"%s",s.c_str());
    fclose(fp);
 
    return true;
}

bool g::io::WriteSignalSequence_int(const char* name, const std::vector<int>& signals)
{
    std::ostringstream o;
    for(size_t i = 0; i < signals.size(); i++){
        o<<signals[i]<<std::endl;
    }
    std::string s=o.str();

    //--- output ----//
    FILE *fp=fopen(name,"wb");
    fprintf(fp,"%s",s.c_str());
    fclose(fp);

    return true;
}


