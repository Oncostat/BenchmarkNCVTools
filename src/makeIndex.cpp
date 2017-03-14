#include <iostream>
#include <string>
#include <fstream>
#include <map>
//#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/unordered_map.hpp>
//#include "vsfunctions.h"
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
    //vector<int> posChrCharSorted{0,11,15,16,17,18,19,20,21,1,2,3,4,5,6,7,8,9,10,12,13,14,22,23,24}; //order of chr in the file if sorted by character, 24 added for the eof
    //vector<int> posChrCharSorted{0,9,10,11,12,13,14,15,16,17,18,1,19,20,21,2,3,4,5,6,7,8,22,23,24};

    //0,9,10,11,12,13,14,15,16,17,18,1,19,20,21,2,3,4,5,6,7,8,22,23;
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



//        cadd.seekg(0, ios::end);dann.seekg(0, ios::end);fath.seekg(0, ios::end);funs.seekg(0, ios::end);
//        streampos caddEof=cadd.tellg(),dannEof=dann.tellg(),fathEof=fath.tellg(),funsEof=funs.tellg(); //Getting the length of file or positon of enf of file
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

        //for(int i=0;i<24;++i){cout << i << "  " << posChrDann.at(i)<< "  vs  " << dannVecChr.at(i) << "  " << posChrDann.at(dannVecChr.at(i)) << endl;}
//        for(int i=0;i<24;++i){cout << posChrDann.at(i) << '\t';}
//        cout << endl;

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
//            cout << i << "  " << saveVecTmp.at(1)<< "  vs  " << dannVecChr.at(i) << "  " << posChrDann.at(dannVecChr.at(i)) << endl;

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
//                cout << "go for " << nextGenPos<<endl;
                vecTmp.at(0)=getFilePos(nextGenPos,cadd,initCadd.zeroBased,posChrCadd.at(caddVecChr.at(i)),initCadd.end,saveVecTmp.at(0));
//                cout << "cadd ok"<<endl;
                //cout << i << "  " << chrSize.size() << "   " << dannVecChr.size() << endl;
/*                if(i==1){
                    //cout << "test chrsize   " << dannVecChr.at(i) << "   " << chrSize.at(dannVecChr.at(i)) <<endl;
                    cout << "test posnextchr   " << dannVecChr.at(i+1) << "   " << posChrDann.at(dannVecChr.at(i+1)) <<endl;
                }*/
                vecTmp.at(1)=getFilePos(nextGenPos,dann,initDann.zeroBased,posChrDann.at(dannVecChr.at(i)),initDann.end,saveVecTmp.at(1));
//                cout << "dann ok"<<endl;
                vecTmp.at(2)=getFilePos(nextGenPos,fath,initFath.zeroBased,posChrFath.at(fathVecChr.at(i)),initFath.end,saveVecTmp.at(2));
//                cout << "fath ok"<<endl;
                vecTmp.at(3)=getFilePos(nextGenPos,funs,initFuns.zeroBased,posChrFuns.at(funsVecChr.at(i)),initFuns.end,saveVecTmp.at(3));
//                cout << "funs ok"<<endl;
                for(int j=0;j<4;++j){
                    if(vecTmp.at(j)==saveVecTmp.at(j)){
                        ++test;
                    }
                }
/*                if(i==1){
                    cout << vecTmp.at(1) << endl;
                    cout << vecTmp.size() << endl;
                    cout << vecTmp.at(vecTmp.size()-1) << endl;
                    cout << vecTmp.at(vecTmp.size()-2) << endl;
                    cout << test << endl;
                    cout<< "test" << endl;
                }*/
                if(test==4 && vecPosTmp.genPos.size()!=0){//same seek position in all files, replace by higher threshold
/*                    if(i==12){
                        cout << "here1  " << nextGenPos<<endl;
                        cout << vecPosTmp.genPos.size() << endl;
                        cout << vecPosTmp.genPos.size()-1<<endl;
                    }*/
                    vecPosTmp.genPos.at(vecPosTmp.genPos.size()-1)=nextGenPos;
                }else{ //add new pos
//                    if(i==12){cout << "here2"<<endl;}
                    vecPosTmp.genPos.push_back(nextGenPos);
                    vecPosTmp.filePos.push_back(vecTmp);
                }
                //if(i==1){cout << "test" <<endl;}

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

/*        for(int i=0;i<24;++i){
            //chrCoord

            while(nextGenPos<chrSize.at(i)){
//use chr coord ?
                posTmp.push_back(indexing(cadd,initCadd,i,nextGenPos,caddEof)); //search filePos of nextGenPos
                posTmp.push_back(indexing(dann,initDann,i,nextGenPos,dannEof));
                posTmp.push_back(indexing(fath,initFath,i,nextGenPos,fathEof));
                posTmp.push_back(indexing(funs,initFuns,i,nextGenPos,funsEof));
                test=0;
                currentSize=vecTmp.filePos.size();
                for(int j=0;j<4;++j){
                    if(posTmp.at(j)==saveInd.at(j)){
                        ++test;
                    }
                }
                if(test!=4){
                    vecTmp.genPos.push_back(nextGenPos);
                    vecTmp.filePos.push_back(posTmp);
                    vecTmp.classified.push_back(fillTmp);
                    savePosTmp=posTmp;
                }
                nextGenPos+=stepGenPos;
            }
            out.push_back(vecTmp);
        }

        for(int i=0;i<out.size();++i){
            for(int j=0;j<out[i].genPos.at(0).size();++j){
                index << i << '\t' << out[i].genPos.at(j);
                for(int k=0;k<out[i].filePos.size();++k){
                    index << '\t' << out[i].filePos.at(j).at(k)<<'\n';
                }
            }
        }
        cout << "Index file was created."<<endl;
    }*/

/*    cout << initCadd.begin << '\t' << initCadd.end << '\t' << initCadd.chrIndic << '\t' << initCadd.zeroBased << endl;
    cout << initDann.begin << '\t' << initDann.end << '\t' << initDann.chrIndic << '\t' << initDann.zeroBased << endl;
    cout << initFath.begin << '\t' << initFath.end << '\t' << initFath.chrIndic << '\t' << initFath.zeroBased << endl;
    cout << initFuns.begin << '\t' << initFuns.end << '\t' << initFuns.chrIndic << '\t' << initFuns.zeroBased << endl;*/

//    cout << cadd.tellg()<<endl;
//    cout << dann.tellg() << endl;
//    cout << initDann.end << endl;
    //dann.seekg(0);
//    dann.seekg(10, ios::end);
//    cout << dann.tellg()<<endl;


/*    fstream caddIndex("Scores/CADD_scores.index", std::ios::in);
    fstream dannIndex("Scores/DANN_scores.index", std::ios::in);
    fstream fathIndex("Scores/FATH_scores.index", std::ios::in);
    fstream funsIndex("Scores/Funseq2_scores.index", std::ios::in);



    vector<vecPos> tmpVec;
    if(!caddIndex){
        fstream caddIndex("Scores/CADD_scores.index", ios::in | ios::out | ios::trunc);
        if(!caddIndex){cerr << "The CADD index file cannot be created." << std::endl;return 1;}
        cout << "CADD index file could not be found."<<endl;
        cout << "Creating CADD index file. Please wait, it could take few minutes."<<endl;
        tmpVec=indexing(cadd,initCadd);
        for(int i=0;i<tmpVec.size();++i){
            for(int j=0;j<tmpVec[i].filePos.size();++j){
                caddIndex << i << '\t' << tmpVec[i].genPos.at(j) << '\t' << tmpVec[i].filePos.at(j)<<endl;
            }
        }
        cout << "CADD index file was created."<<endl;
    }
    if(!dannIndex){
        fstream dannIndex("Scores/DANN_scores.index", ios::in | ios::out | ios::trunc);
        if(!dannIndex){cerr << "The DANN index file cannot be created." << std::endl;return 1;}
        cout << "DANN index file could not be found."<<endl;
        cout << "Creating DANN index file. Please wait, it could take few minutes."<<endl;
        tmpVec=indexing(dann,initDann);
        for(int i=0;i<tmpVec.size();++i){
            for(int j=0;j<tmpVec[i].filePos.size();++j){
                dannIndex << i << '\t' << tmpVec[i].genPos.at(j) << '\t' << tmpVec[i].filePos.at(j)<<endl;
            }
        }
        cout << "DANN index file was created."<<endl;
    }
    if(!fathIndex){
        fstream fathIndex("Scores/FATH_scores.index", ios::in | ios::out | ios::trunc);
        if(!fathIndex){cerr << "The FATH index file cannot be created." << std::endl;return 1;}
        cout << "FATHMM-MKL index file could not be found."<<endl;
        cout << "Creating FATHMM-MKL index file. Please wait, it could take few minutes."<<endl;
        tmpVec=indexing(fath,initFath);
        for(int i=0;i<tmpVec.size();++i){
            for(int j=0;j<tmpVec[i].filePos.size();++j){
                fathIndex << i << '\t' << tmpVec[i].genPos.at(j) << '\t' << tmpVec[i].filePos.at(j)<<endl;
            }
        }
        cout << "FATHMM-MKL index file was created."<<endl;
    }
    if(!funsIndex){
        fstream funsIndex("Scores/Funseq2_scores.index", ios::in | ios::out | ios::trunc);
        if(!funsIndex){cerr << "The Funseq2 index file cannot be created." << std::endl;return 1;}
        cout << "Funseq2 index file could not be found."<<endl;
        cout << "Creating Funseq2 index file. Please wait, it could take few minutes."<<endl;
        tmpVec=indexing(funs,initFuns);
        for(int i=0;i<tmpVec.size();++i){
            for(int j=0;j<tmpVec[i].filePos.size();++j){
                funsIndex << i << '\t' << tmpVec[i].genPos.at(j) << '\t' << tmpVec[i].filePos.at(j)<<endl;
            }
        }
        cout << "Funseq2 index files was created."<<endl;
    }
*/
