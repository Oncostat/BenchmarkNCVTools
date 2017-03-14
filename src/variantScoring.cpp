#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/unordered_map.hpp>
#include "vsfunctions.h"

using namespace std;




int main(int argc, char *argv[]){
    int test=0;

    fstream data(argv[1], std::ios::in);
    if(!data){cerr << "The data file does not exist." << endl;return 1;}
    ofstream output(argv[2], ios::out | ios::trunc);
    if(!output){cerr << "The output file cannot be created." << std::endl;return 1;}

    fstream cadd("Scores/CADD_scores.tsv", std::ios::in);
    fstream dann("Scores/DANN_scores.tsv", std::ios::in);
    fstream fath("Scores/FATH_scores.tsv", std::ios::in);
    fstream funs("Scores/Funseq2_scores.tsv", std::ios::in);

    if(!cadd){cerr << "The CADD score file does not exist." << endl;return 1;}
    if(!dann){cerr << "The DANN score file does not exist." << endl;return 1;}
    if(!fath){cerr << "The FATHMM-MKL score file does not exist." << endl;return 1;}
    if(!funs){cerr << "The FUNSEQ 2 score file does not exist." << endl;return 1;}

    fstream scoresIndex("Scores/index", std::ios::in);
    if(!scoresIndex){cerr << "The index file does not exist." << endl;return 1;}


    fstream initFile("Scores/initScores", std::ios::in);
    if(!initFile){cerr << "The score file for algorithm initiation does not exist." << endl;return 1;}


    struct fileInfos initData=initializeData(data,"tsv");
    load(data,(streampos)0);

    struct fileInfos initCadd=fromInit(initFile);
    struct fileInfos initDann=fromInit(initFile);
    struct fileInfos initFath=fromInit(initFile);
    struct fileInfos initFuns=fromInit(initFile);

    vector<vecPos> index=loadIndex(scoresIndex);


    load(data,initData.begin);
    load(cadd,initCadd.begin);
    load(dann,initDann.begin);
    load(fath,initFath.begin);
    load(funs,initFuns.begin);


    string lData;
    struct position posRef;
    vector<string> file;
    int countIn=0,saveChr=0;
    long long int saveGenPos=saveGenPos=index.at(0).genPos.at(0),countGenPos,savePos=0;
    while(getline(data,lData)){
        file.push_back(lData);     //Data loading
        ++countIn;
        if(countIn==5000000){
            for(int i=0;i<file.size();++i){
                posRef=posStruct(file.at(i),initData.chrIndic,initData.zeroBased);
                if(posRef.chr!=saveChr || posRef.posDNA>saveGenPos){
                    if(posRef.chr!=saveChr){
                        if(posRef.chr>=index.size()){
                            file.at(i)+="\tNA\tNA\tNA\tNA\n";
                            continue;
                        }
                        saveChr=posRef.chr;
                        countGenPos=0;
                        saveGenPos=index.at(saveChr).genPos.at(0);
                    }
                    while(posRef.posDNA>saveGenPos){
                        ++countGenPos;
                        saveGenPos=index.at(saveChr).genPos.at(countGenPos);
                    }
                    if(countGenPos==0){
                        file.at(i)+="\tNA\tNA\tNA\tNA\n";
                        continue;
                    }else{
                        load(cadd,index.at(saveChr).filePos.at(countGenPos-1).at(0));
                        load(dann,index.at(saveChr).filePos.at(countGenPos-1).at(1));
                        load(fath,index.at(saveChr).filePos.at(countGenPos-1).at(2));
                        load(funs,index.at(saveChr).filePos.at(countGenPos-1).at(3));
                    }
                }
                scoring(&posRef,cadd,initCadd.chrIndic,initCadd.zeroBased,&file.at(i));
                scoring(&posRef,dann,initDann.chrIndic,initDann.zeroBased,&file.at(i));
                scoring(&posRef,fath,initFath.chrIndic,initFath.zeroBased,&file.at(i));
                scoring(&posRef,funs,initFuns.chrIndic,initFuns.zeroBased,&file.at(i));
                file.at(i)+='\n';
                savePos=posRef.posDNA;
            }
            countIn=0;
            for(int i=0;i<file.size();++i){
                output << file.at(i) << '\n';
            }
            file.clear();
       }
    }


    for(int i=0;i<file.size();++i){
        posRef=posStruct(file.at(i),initData.chrIndic,initData.zeroBased);
        if(posRef.chr!=saveChr || posRef.posDNA>saveGenPos){
            if(posRef.chr!=saveChr){
                if(posRef.chr>=index.size()){
                    file.at(i)+="\tNA\tNA\tNA\tNA\n";
                    continue;
                }
                saveChr=posRef.chr;
                countGenPos=0;
                saveGenPos=index.at(saveChr).genPos.at(0);
            }
            while(posRef.posDNA>saveGenPos){
                ++countGenPos;
                saveGenPos=index.at(saveChr).genPos.at(countGenPos);
            }
            if(countGenPos==0){
                file.at(i)+="\tNA\tNA\tNA\tNA\n";
                continue;
            }else{
                load(cadd,index.at(saveChr).filePos.at(countGenPos-1).at(0));
                load(dann,index.at(saveChr).filePos.at(countGenPos-1).at(1));
                load(fath,index.at(saveChr).filePos.at(countGenPos-1).at(2));
                load(funs,index.at(saveChr).filePos.at(countGenPos-1).at(3));
            }
        }
        scoring(&posRef,cadd,initCadd.chrIndic,initCadd.zeroBased,&file.at(i));
        scoring(&posRef,dann,initDann.chrIndic,initDann.zeroBased,&file.at(i));
        scoring(&posRef,fath,initFath.chrIndic,initFath.zeroBased,&file.at(i));
        scoring(&posRef,funs,initFuns.chrIndic,initFuns.zeroBased,&file.at(i));
        file.at(i)+='\n';
    }


    countIn=0;
    for(int i=0;i<file.size();++i){
        output << file.at(i);
    }
    file.clear();


    data.close();
    cadd.close();
    dann.close();
    fath.close();
    funs.close();
    output.close();

    return 0;
}
