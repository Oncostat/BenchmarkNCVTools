#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <boost/algorithm/string.hpp>



using namespace std;

struct rec{
    vector<bool> recConf;
    vector<bool> recUnkn;
    vector<bool> recTotal;
};



void seek(ifstream &file, streampos pos){
    file.clear();
    file.seekg(pos);
}

int main(int argc, char* argv[]){
	int minrec,maxrec;
        ifstream file(argv[1], std::ios::in);
        if(!file){cerr << "The input file does not exist." << endl;return 1;}

	if(argc==2){
		minrec=0;
		maxrec=20;
		cout << "Warnings ! No recurrence thresholds were provided. The thresholds were taken from 1 to 20." << endl;
	}else{
		if(argc!=4){
			cout << "Please indicate the minimum  and maximum recurrence thresholds (or let both unspecified for range from 1 to 20)." << endl;
			return 1;
		}else{
                        minrec=stoi(argv[2]);
                        maxrec=stoi(argv[3]);
		}
	}
	if(minrec>maxrec){
                cout << "Program was stopped because minimum recurrence threshold > maximum recurrence threshold.\nPlease correctly provide the arguments (inputFile minimalThreshold maximalThreshold)" << endl;
		return 1;
	}

        vector<int> recs;
        for(int i=minrec;i<=maxrec;++i){recs.push_back(i);}
        string lfile;
        int recConf,recUnkn,recTotal;


        rec out;
        vector<rec> outs;
        vector<string> data;
        while(getline(file,lfile,'\n')){data.push_back(lfile);}
        file.close();
        int maxConfRec=0,maxUnknRec=0,maxTotalRec=0;

        for(int i=minrec;i<=maxrec;++i){
            for(int j=0;j<data.size();++j){
                vector<string> strs;
                boost::split(strs, data.at(j), boost::is_any_of("\t"));
                recConf = stoi(strs.at(5));
                recUnkn = stoi(strs.at(6));
                recTotal = stoi(strs.at(7));
                if(recConf>=i){maxConfRec=i;out.recConf.push_back(true);}else{out.recConf.push_back(false);}
                if(recUnkn>=i){maxUnknRec=i;out.recUnkn.push_back(true);}else{out.recUnkn.push_back(false);}
                if(recTotal>=i){maxTotalRec=i;out.recTotal.push_back(true);}else{out.recTotal.push_back(false);}
            }
            outs.push_back(out);
            out.recConf.clear();
            out.recUnkn.clear();
            out.recTotal.clear();
        }
        for(int i=0;i<recs.size();++i){
            ofstream outFileConf("Databases/COSMICrecConf"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
            if(!outFileConf){cerr << "The output files Databases/COSMICrecConf" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
            for(int j=0;j<data.size();++j){
                if(outs.at(i).recConf.at(j)){outFileConf << data.at(j) << '\n';}
            }
            outFileConf.close();
            ofstream outFileUnkn("Databases/COSMICrecUnkn"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
            if(!outFileUnkn){cerr << "The output files Databases/COSMICrecUnkn" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
            for(int j=0;j<data.size();++j){
                if(outs.at(i).recUnkn.at(j)){outFileUnkn << data.at(j) << '\n';}
            }
            outFileUnkn.close();
            ofstream outFileTotal("Databases/COSMICrecTotal"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
            if(!outFileTotal){cerr << "The output files Databases/COSMICrecTotal" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
            for(int j=0;j<data.size();++j){
                if(outs.at(i).recTotal.at(j)){outFileTotal << data.at(j) << '\n';}
            }
            outFileTotal.close();
        }



	return 0;
}

