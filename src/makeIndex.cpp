#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/unordered_map.hpp>
#include "makeIndex.h"


using namespace std;


int main(){
    fstream cadd("Scores/CADD_scores.tsv", std::ios::in);
    fstream dann("Scores/DANN_scores.tsv", std::ios::in);
    fstream fath("Scores/FATH_scores.tsv", std::ios::in);
    fstream funs("Scores/Funseq2_scores.tsv", std::ios::in);

    if(!cadd){cerr << "The CADD score file does not exist." << endl;return 1;}
    if(!dann){cerr << "The DANN score file does not exist." << endl;return 1;}
    if(!fath){cerr << "The FATHMM-MKL score file does not exist." << endl;return 1;}
    if(!funs){cerr << "The FUNSEQ 2 score file does not exist." << endl;return 1;}


    struct fileInfos initCadd=initializeData(cadd,"Scores/CADD_scores.tsv");
    struct fileInfos initDann=initializeData(dann,"Scores/DANN_scores.tsv");
    struct fileInfos initFath=initializeData(fath,"Scores/FATH_scores.tsv");
    struct fileInfos initFuns=initializeData(funs,"Scores/Funseq2_scores.tsv");

    vector<int> posChr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    vector<int> posChrCharSorted{9,19,3,4,5,6,7,8,22,10,11,12,13,14,15,16,17,18,1,20,21,2,23,24};

    vector<int> caddVecChr,dannVecChr,fathVecChr,funsVecChr;
    if(initCadd.charSorted){caddVecChr=posChrCharSorted;}else{caddVecChr=posChr;}
    if(initDann.charSorted){dannVecChr=posChrCharSorted;}else{dannVecChr=posChr;}
    if(initFath.charSorted){fathVecChr=posChrCharSorted;}else{fathVecChr=posChr;}
    if(initFuns.charSorted){funsVecChr=posChrCharSorted;}else{funsVecChr=posChr;}

    fstream initFile("Scores/initScores", std::ios::in);
    if(!initFile){
        ofstream initFile("Scores/initScores", ios::out | ios::trunc);
        if(!initFile){cerr << "The score file for algorithm initiation cannot be created." << std::endl;return 1;}
        initFile << "CADD\t" << initCadd.begin << '\t' << initCadd.end << '\t' << initCadd.chrIndic << '\t' << initCadd.zeroBased << "\n";
        initFile << "DANN\t" << initDann.begin << '\t' << initDann.end << '\t' << initDann.chrIndic << '\t' << initDann.zeroBased << "\n";
        initFile << "FATH\t" << initFath.begin << '\t' << initFath.end << '\t' << initFath.chrIndic << '\t' << initFath.zeroBased << "\n";
        initFile << "Funs\t" << initFuns.begin << '\t' << initFuns.end << '\t' << initFuns.chrIndic << '\t' << initFuns.zeroBased << "\n";

    }


    fstream index("Scores/index", std::ios::in);
    if(!index){
        ofstream index("Scores/index", ios::out | ios::trunc);
        if(!index){cerr << "The index file cannot be created." << std::endl;return 1;}
        cout << "Index file could not be found."<<endl;
        cout << "Creating index file. Please wait, it could take few minutes."<<endl;



        struct vecPos vecPosTmp;
        vector<vecPos> out(24);
        vector<streampos> posChrCadd,posChrDann,posChrFath,posChrFuns;
        vector<long long int> chrSize{249250622,243199374,198022431,191154277,180915261,171115068,159138664,146364023,141213432,135534748,135006517,133851896,115169879,107349541,102531393,90354754,81195211,78077249,59128984,63025521,48129896,51304567,155270561,59373567};
        long long int nextGenPos=1;


        posChrCadd=indexChr(cadd,initCadd,initCadd.end);
        posChrDann=indexChr(dann,initDann,initDann.end);
        posChrFath=indexChr(fath,initFath,initFath.end);
        posChrFuns=indexChr(funs,initFuns,initFuns.end);

        vector<streampos> vecTmp(4),saveVecTmp(4);
        int test;
        int chrDone,progress=0;


        for(int i=0;i<24;++i){
            cout << "chr" << i+1<<endl;
            cout.flush();
            if(initCadd.charSorted){load(cadd,posChrCadd.at(i));}
            if(initDann.charSorted){load(dann,posChrDann.at(i));}
            if(initFath.charSorted){load(fath,posChrFath.at(i));}
            if(initFuns.charSorted){load(funs,posChrFuns.at(i));}
            saveVecTmp.at(0)=posChrCadd.at(i);
            saveVecTmp.at(1)=posChrDann.at(i);
            saveVecTmp.at(2)=posChrFath.at(i);
            saveVecTmp.at(3)=posChrFuns.at(i);

            nextGenPos=1;
            while(nextGenPos<chrSize.at(i)){
                test=0;
                nextGenPos+=100000;
                chrDone=(int)round(100*nextGenPos/chrSize.at(i));
                if(chrDone!=progress && (chrDone%5)==0){
                    cout << chrDone << "%\t";
                    cout.flush();
                    progress=chrDone;
                }
                vecTmp.at(0)=getFilePos(nextGenPos,cadd,initCadd.zeroBased,posChrCadd.at(caddVecChr.at(i)),initCadd.end,saveVecTmp.at(0));
                vecTmp.at(1)=getFilePos(nextGenPos,dann,initDann.zeroBased,posChrDann.at(dannVecChr.at(i)),initDann.end,saveVecTmp.at(1));
                vecTmp.at(2)=getFilePos(nextGenPos,fath,initFath.zeroBased,posChrFath.at(fathVecChr.at(i)),initFath.end,saveVecTmp.at(2));
                vecTmp.at(3)=getFilePos(nextGenPos,funs,initFuns.zeroBased,posChrFuns.at(funsVecChr.at(i)),initFuns.end,saveVecTmp.at(3));
                for(int j=0;j<4;++j){
                    if(vecTmp.at(j)==saveVecTmp.at(j)){
                        ++test;
                    }
                }
                if(test==4 && vecPosTmp.genPos.size()!=0){//same seek position in all files, replace by higher threshold
                    vecPosTmp.genPos.at(vecPosTmp.genPos.size()-1)=nextGenPos;
                }else{ //add new pos
                    vecPosTmp.genPos.push_back(nextGenPos);
                    vecPosTmp.filePos.push_back(vecTmp);
                }
                saveVecTmp=vecTmp;
            }
            out.at(i)=vecPosTmp;
            vecPosTmp.genPos.clear();
            vecPosTmp.filePos.clear();
            cout << endl;
        }
        for(int i=0;i<out.size();++i){
            for(int j=0;j<out[i].genPos.size();++j){
                index << i << '\t' << out[i].genPos.at(j) << '\t';
                for(int k=0;k<4;++k){
                    index << out[i].filePos.at(j).at(k) << '\t';
                }
                index << '\n';
            }
        }
        cout << "Index file was created."<<endl;
    }

    cadd.close();
    dann.close();
    fath.close();
    funs.close();

    return 0;
}


