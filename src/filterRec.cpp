/*
 * filter_rec.cpp
 *
 *  Created on: 18 déc. 2015
 *      Author: drubay
 */



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
//    string line;
//    vector<bool> recConf,recUnkn,recTotal;
    vector<bool> recConf;
    vector<bool> recUnkn;
    vector<bool> recTotal;
};



void seek(ifstream &file, streampos pos){
    file.clear();
    file.seekg(pos);
}

int main(int argc, char* argv[]){
	int minrec,maxrec; //,reccheck
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
        int recConf,recUnkn,recTotal;//,nbRec=minrec-maxrec+1


        rec out;
        vector<rec> outs;
        //if(argv[1].find("conf")!=-1){int labels=1;}else{int labels=0;}
        vector<string> data;
        while(getline(file,lfile,'\n')){data.push_back(lfile);}
        file.close();
        int maxConfRec=0,maxUnknRec=0,maxTotalRec=0;

        //search the max rec of each categ + vector saving the number of the lines with rec
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


        /* for each line
         *
         */


/*            vector<string> strs;
            boost::split(strs, ldata, boost::is_any_of("\t"));
            recConf = stoi(strs.at(5));
            recUnkn = stoi(strs.at(6));
            recTotal = stoi(strs.at(7));
            out.line=ldata;
            for(int i=minrec;i<=maxrec;++i){
                if(recConf>=i){out.recConf.push_back(true);++countRecConf.at(i-minrec);}else{out.recConf.push_back(false);}
                if(recUnkn>=i){out.recUnkn.push_back(true);++countRecUnkn.at(i-minrec);}else{out.recUnkn.push_back(false);}
                if(recTotal>=i){out.recTotal.push_back(true);++countRecTotal.at(i-minrec);}else{out.recTotal.push_back(false);}
            }
            outs.push_back(out);
        }
        int numRec;
        for(int i=0;i<nbRec;++i){
            numRec=recs.at(i);
            if(countRecConf.at(i)>0){
                ofstream outConf("Databases/COSMICrecConf"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
                if(!outConf){cerr << "The output files Databases/COSMICrecConf" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
            }
            if(countRecUnkn.at(i)>0){
                ofstream outUnkn("Databases/COSMICrecUnkn"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
                if(!outUnkn){cerr << "The output files Databases/COSMICrecUnkn" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
            }
            if(countRecTotal.at(i)>0){
                ofstream outTotal("Databases/COSMICrecTotal"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
                if(!outTotal){cerr << "The output files Databases/COSMICrecTotal" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
            }
            for(int j=0;j<outs.size();++j){
                if(outs.at(j).recConf.at(i)){

                    for(int j=0;j<outs.at(i).outConf.size();++j){seek(data,outs.at(i).outConf.at(j));getline(data,ldata);outFile << outs.at(i).outConf.at(j) << endl;}
                    outFile.close();
                }
            }
            if(countRecConf.at(i)>0){outUnkn.close();}
            if(countRecUnkn.at(i)>0){outTotal.close();}
            if(countRecTotal.at(i)>0){outConf.close();}


        }*/

	return 0;
}

/*            ofstream conf("Databases/COSMICrecConf"+to_string(i)+".txt",ios::out | ios::trunc);
            ofstream unkn("Databases/COSMICrecUnkn"+to_string(i)+".txt",ios::out | ios::trunc);
            ofstream total("Databases/COSMICrecTotal"+to_string(i)+".txt",ios::out | ios::trunc);
            if(!conf){cerr << "The output files cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}*/
/*                if(labels==0){
                        tmp = "Databases/COSMICrec" + to_string(i) + ".txt";
                }else{
                        tmp = "Databases/COSMICconfrec" + to_string(i) + ".txt";
                }*/

//                outputFiles.push_back(tmp);
                //targv[1]est output and remove, by the way, old files

//                                tmpData.open(outputFilesConf[i-minrec], ios_base::app | ios::out);
//                                tmpData << ldata << endl;
//                                tmpData.close();
//                                tmpData.open(outputFilesUnkn[i-minrec], ios_base::app | ios::out);
//                                tmpData << ldata << endl;
//                                tmpData.close();
//                                tmpData.open(outputFilesTotal[i-minrec], ios_base::app | ios::out);
//                                tmpData << ldata << endl;
//                                tmpData.close();




/*
if(maxConfRec>=recs.at(i)){
    ofstream outFile("Databases/COSMICrecConf"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
    if(!outFile){cerr << "The output files Databases/COSMICrecConf" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
    for(int j=0;j<data.size();++j){
        if(outs.at(i).recConf.at(j)){outFile << data.at(j) << '\n';}
    }
    outFile.close();
}
if(maxUnknRec>=recs.at(i)){
    ofstream outFile("Databases/COSMICrecUnkn"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
    if(!outFile){cerr << "The output files Databases/COSMICrecUnkn" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
    for(int j=0;j<data.size();++j){
        if(outs.at(i).recUnkn.at(j)){outFile << data.at(j) << '\n';}
    }
    outFile.close();
}
if(maxTotalRec>=recs.at(i)){
    ofstream outFile("Databases/COSMICrecTotal"+to_string(recs.at(i))+".txt",ios::out | ios::trunc);
    if(!outFile){cerr << "The output files Databases/COSMICrecTotal" << to_string(recs.at(i)) << ".txt cannot be created.\nPlease check that the \"Databases\" folder is in the same folder than the snakefile." << endl;return 1;}
    for(int j=0;j<data.size();++j){
        if(outs.at(i).recTotal.at(j)){outFile << data.at(j) << '\n';}
    }
    outFile.close();
}
*/
