
#include <iostream>
#include <fstream>
#include <istream>
#include <stdio.h>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;


int main(int argc, char* argv[]){
    ifstream file(argv[1], std::ios::in);
    if(!file){cerr << "The input file does not exist." << std::endl;return 1;}
    ofstream ouput(argv[2], ios::out | ios::trunc);
    if(!ouput){cerr << "The output file cannot be created." << std::endl;return 1;}
    ifstream features("Databases/feat_gencode.txt", std::ios::in);
    if(!features){cerr << "The feature file does not exist." << endl;return 1;}

    string lout;
    vector<string> out;
    lout="chr\tpos";
    string lfile;
    vector<string> feats;
    while(getline(features,lfile)){feats.push_back(lfile);lout+='\t'+lfile;}
    feats.push_back("-1");lout+="\t-1";
    out.push_back(lout);

    vector<string> genPos(2);
    vector<int> featureCount(feats.size(),0);
    vector<string> strs;
    vector<string>::iterator itFeat;

    getline(file,lfile);
    boost::split(strs, lfile, boost::is_any_of("\t"));
    genPos[0]=strs[0];
    genPos[1]=strs[1];
    for(int i=2;i<strs.size();++i){
        strs[i].erase(remove_if(strs[i].begin(),strs[i].end(),::isspace),strs[i].end());
        itFeat=find(feats.begin(),feats.end(),strs[i]);
        if(itFeat!=feats.end()){
            auto pos=itFeat-feats.begin();
            ++featureCount[pos];
        }else{

        }
    }

    while(getline(file,lfile)){
        boost::split(strs, lfile, boost::is_any_of("\t"));
        if(genPos[0]!=strs[0] || genPos[1]!=strs[1]){
            lout=genPos[0]+'\t'+genPos[1];
            for(int i=0;i<featureCount.size();++i){
                lout+='\t'+to_string(featureCount[i]); 
            }
            out.push_back(lout);
            fill(featureCount.begin(),featureCount.end(),0); 
            genPos[0]=strs[0];
            genPos[1]=strs[1];
        }
        for(int i=2;i<strs.size();++i){
            strs[i].erase(remove_if(strs[i].begin(),strs[i].end(),::isspace),strs[i].end());
            itFeat=find(feats.begin(),feats.end(),strs[i]);
            if(itFeat!=feats.end()){
                auto pos=itFeat-feats.begin();
                ++featureCount[pos];
            }
        }
    }
    lout=genPos[0]+'\t'+genPos[1];
    for(int i=0;i<featureCount.size();++i){
        lout+='\t'+to_string(featureCount[i]);
    }
    out.push_back(lout);
    for(int i=0;i<out.size();++i){
        ouput << out[i]<<'\n';
    }


    return 0;
}
