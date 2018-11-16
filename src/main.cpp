#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include "evaluate.h"
#include "decimation.h"

EvaluateInfo parseLine(const std::string &line)
{
    EvaluateInfo ret;
    std::string token;
    std::stringstream ss(line);

    ss >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if(token == "glc") {
        ret.algorithm = EvaluateInfo::GLC;
    } else if(token == "none") {
        ret.algorithm = EvaluateInfo::None;
    } else {
        ret.algorithm = EvaluateInfo::NFR;
    }

    ss >> token;
    ret.g2oname = token;

    ss >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if(token == "online") {
        ret.decimate = onlineDecimate;
    } else if(token == "cluster") {
        ret.decimate = clusterDecimate;
    } else /* if(token == "global") */ {
        ret.decimate = globalDecimate;
    }

    ss >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if(token == "tree") {
        ret.sparsityOptions.topology = SparsityOptions::Tree;
    } else if(token == "subgr") {
        ret.sparsityOptions.topology = SparsityOptions::Subgraph;
    } else if(token == "clsubgr") {
        ret.sparsityOptions.topology = SparsityOptions::CliqueySubgraph;
    } else if(token == "dense") {
        ret.sparsityOptions.topology = SparsityOptions::Dense;
    } else /* if(token == "cldense") */ {
        ret.sparsityOptions.topology = SparsityOptions::CliqueyDense;
    }

    ss >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if(token == "local") {
        ret.sparsityOptions.linPoint = SparsityOptions::Local;
    } else /* if(token == "global") */ {
        ret.sparsityOptions.linPoint = SparsityOptions::Global;
    }

    ss >> ret.decimateOptions.sparsity;
    if(ss.good()) {
        ss >> ret.kldPeriod;
    } else {
        ret.kldPeriod = 10;
    }

    if(ret.decimate == globalDecimate) {
        ret.kldPeriod = std::numeric_limits<int>::max();
    }

    if(ss.good()) {
        ss >> token;
        std::transform(token.begin(), token.end(), token.begin(), ::tolower);
        if(token == "chi2") {
            ret.useChi2 = true;
        } else {
            ret.useChi2 = false;
        }
    } else {
        ret.useChi2 = false;
    }

    if(ss.good()) {
        ss >> ret.decimateOptions.clusterSize;
    } else {
        ret.decimateOptions.clusterSize = 100;
    }

    return ret;
}

std::string trim(const std::string &s)
{
    int i = 0, j = s.length();
    while(i < s.length() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n')) i++;
    while(j > 0 && (s[j-1] == ' ' || s[j-1] == '\t' || s[j-1] == '\n')) j--;
    if(i < j) {
        return s.substr(i, j - i);
    } else {
        return std::string();
    }
}

std::vector<EvaluateInfo> loadEvaluateInfo(const char *destdir, const char *infoname)
{
    std::vector<EvaluateInfo> ret;
    std::string line;
    std::ifstream f(infoname);

    while(f.good()) {
        std::getline(f, line);
        line = trim(line);
        if(line.length() == 0 || line[0] == '#') continue;
        ret.push_back(parseLine(line));
        ret.back().destdir = destdir;
    }

    f.close();

    return ret;
}

int main(int argc, char *argv[])
{
    if(argc >= 3) {
        int threads = 0;
        if(argc >= 4) threads = atoi(argv[3]);
        parallelEvaluate(loadEvaluateInfo(argv[1], argv[2]), threads);
    }
    return 0;
}
