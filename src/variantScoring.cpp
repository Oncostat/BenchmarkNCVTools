#include <iostream>
#include <string>
#include <fstream>
#include <map>
//#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/unordered_map.hpp>
#include "vsfunctions.h"
//#include "makeindex.h"

using namespace std;


/*Buffer version
 * HDD block size: 4096
 * control input
 * load segments delimited by posLim when required
 *
 *
 */



int main(int argc, char *argv[]){
    int test=0;
//    cout << "test" <<test++ << endl;

    fstream data(argv[1], std::ios::in);
    if(!data){cerr << "The data file does not exist." << endl;return 1;}
    //ifstream data("Databases/Clinvar.gw", std::ios::in);
//    ifstream data("Databases/COSMIC.gw", std::ios::in);
    //ofstream output("output.tsv", ios::out | ios::trunc);
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
//                    cout << "change  " << file.at(i)<<endl;
                    if(posRef.chr!=saveChr){
                        if(posRef.chr>=index.size()){
                            file.at(i)+="\tNA\tNA\tNA\tNA\n";
                            continue;
                        }
                        saveChr=posRef.chr;
                        countGenPos=0;
                        saveGenPos=index.at(saveChr).genPos.at(0);
                    }
                    while(posRef.posDNA>saveGenPos){ //find closest (higher) pos in index
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
        //cout << file.at(i)<< endl;
        if(posRef.chr!=saveChr || posRef.posDNA>saveGenPos){
            if(posRef.chr!=saveChr){
                if(posRef.chr>=index.size()){
                    file.at(i)+="\tNA\tNA\tNA\tNA\n";
                    continue;
                }
                saveChr=posRef.chr;
                countGenPos=0;
                saveGenPos=index.at(saveChr).genPos.at(0);
//                cout << saveChr << '\t'<<saveGenPos<<endl;
            }
            while(posRef.posDNA>saveGenPos){ //find closest (higher) pos in index
                ++countGenPos;
                saveGenPos=index.at(saveChr).genPos.at(countGenPos);
            }
            if(countGenPos==0){
                file.at(i)+="\tNA\tNA\tNA\tNA\n";
                continue;
            }else{
//                cout << saveChr << '\t' << countGenPos-1 << endl;
//                cout << "genpos:  " << index.at(saveChr).genPos.at(countGenPos-1) << endl;
//                for(int j=0;j<4;++j){cout << index.at(saveChr).filePos.at(countGenPos-1).at(j)<<'\t';}
                load(cadd,index.at(saveChr).filePos.at(countGenPos-1).at(0));
                load(dann,index.at(saveChr).filePos.at(countGenPos-1).at(1));
                load(fath,index.at(saveChr).filePos.at(countGenPos-1).at(2));
                load(funs,index.at(saveChr).filePos.at(countGenPos-1).at(3));
            }
        }
/*        streampos temp=dann.tellg();
        string line;
        getline(dann,line);
        cout << line << endl;
        dann.clear();
        dann.seekg(temp);*/
//        cout << "test0:  " << file.at(i) <<endl;
        scoring(&posRef,cadd,initCadd.chrIndic,initCadd.zeroBased,&file.at(i));
//        cout << "test1"<<endl;
//        cout << file.at(i) << endl;
        scoring(&posRef,dann,initDann.chrIndic,initDann.zeroBased,&file.at(i));
//        cout << file.at(i) << endl<<endl;
//        cout << "test2"<<endl;
        scoring(&posRef,fath,initFath.chrIndic,initFath.zeroBased,&file.at(i));
//        cout << "test3"<<endl;
        scoring(&posRef,funs,initFuns.chrIndic,initFuns.zeroBased,&file.at(i));
//        cout << "test4"<<endl;
        file.at(i)+='\n';
//        cout << "test5"<<endl;
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


/*    ////////////////////////////////////////////////////////////////////////////
    //CADD
    for(int i=0;i<file.size();++i){
        posRef=posStruct(lData,initData.chrIndic,initData.zeroBased);
        if(posRef.chr!=saveChr){
            saveChr=posRef.chr;
            saveCaddPos=0;
            updateThres(&caddThres,&caddPos,saveChr,saveCaddPos);
        }
        if(posRef.posDNA>caddThres.nextPos){
            while(posRef.posDNA>caddPos.at(saveChr).genPos.at(saveCaddPos+1)){++saveCaddPos;}
            updateThres(&caddThres,&caddPos,saveChr,saveCaddPos);
            load(cadd,caddThres.start);
        }
        pos=cadd.tellg();
        while(getline(cadd,line)){
            posTmp=posStruct(line,initCadd.chrIndic,initCadd.zeroBased);
            if(posTmp.chr==posRef.chr){
                if(posTmp.posDNA==posRef.posDNA){
                    if(posTmp.alt.find(posRef.alt)!=std::string::npos){
                        file.at(i)+='\t'+posTmp.score;
                        break;
                    }else{
                        if(posTmp.alt.size()==1 && posTmp.alt<posRef.alt){
                            continue;
                        }
                    }
                }else{
                    if(posTmp.posDNA<posRef.posDNA){
                        continue;
                    }
                }
            }
            load(cadd,pos);
            file.at(i)+="\tNA";
            break;
        }
    }


    ////////////////////////////////////////////////////////////////////////////
    //DANN
    for(int i=0;i<file.size();++i){
        posRef=posStruct(lData,initData.chrIndic,initData.zeroBased);
        if(posRef.chr!=saveChr){
            saveChr=posRef.chr;
            saveDannPos=0;
            updateThres(&dannThres,&dannPos,saveChr,saveDannPos);
        }
        if(posRef.posDNA>dannThres.nextPos){
            while(posRef.posDNA>dannPos.at(saveChr).genPos.at(saveDannPos+1)){++saveDannPos;}
            updateThres(&dannThres,&dannPos,saveChr,saveDannPos);
            load(dann,dannThres.start);
        }
        pos=dann.tellg();
        while(getline(dann,line)){
            posTmp=posStruct(line,initDann.chrIndic,initDann.zeroBased);
            if(posTmp.chr==posRef.chr){
                if(posTmp.posDNA==posRef.posDNA){
                    if(posTmp.alt.find(posRef.alt)!=std::string::npos){
                        file.at(i)+='\t'+posTmp.score;
                        break;
                    }else{
                        if(posTmp.alt.size()==1 && posTmp.alt<posRef.alt){
                            continue;
                        }
                    }
                }else{
                    if(posTmp.posDNA<posRef.posDNA){
                        continue;
                    }
                }
            }
            load(dann,pos);
            file.at(i)+="\tNA";
            break;
        }
    }


    ////////////////////////////////////////////////////////////////////////////
    //FATH
    for(int i=0;i<file.size();++i){
        posRef=posStruct(lData,initData.chrIndic,initData.zeroBased);
        if(posRef.chr!=saveChr){
            saveChr=posRef.chr;
            saveFathPos=0;
            updateThres(&fathThres,&fathPos,saveChr,saveFathPos);
        }
        if(posRef.posDNA>fathThres.nextPos){
            while(posRef.posDNA>fathPos.at(saveChr).genPos.at(saveFathPos+1)){
                ++saveFathPos;
            }
            updateThres(&fathThres,&fathPos,saveChr,saveFathPos);
            load(fath,fathThres.start);
        }
        pos=fath.tellg();
        while(getline(fath,line)){
            posTmp=posStruct(line,initFath.chrIndic,initFath.zeroBased);
            if(posTmp.chr==posRef.chr){
                if(posTmp.posDNA==posRef.posDNA){
                    if(posTmp.alt.find(posRef.alt)!=std::string::npos){
                        file.at(i)+='\t'+posTmp.score;
                        break;
                    }else{
                        if(posTmp.alt.size()==1 && posTmp.alt<posRef.alt){
                            continue;
                        }
                    }
                }else{
                    if(posTmp.posDNA<posRef.posDNA){
                        continue;
                    }
                }
            }
            load(fath,pos);
            file.at(i)+="\tNA";
            break;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    //Funseq2
    for(int i=0;i<file.size();++i){
        posRef=posStruct(lData,initData.chrIndic,initData.zeroBased);
        if(posRef.chr!=saveChr){
            saveChr=posRef.chr;
            saveFunsPos=0;
            updateThres(&funsThres,&funsPos,saveChr,saveFunsPos);
        }
        if(posRef.posDNA>funsThres.nextPos){
            while(posRef.posDNA>funsPos.at(saveChr).genPos.at(saveFunsPos+1)){++saveFunsPos;}
            updateThres(&funsThres,&funsPos,saveChr,saveFunsPos);
            load(funs,funsThres.start);
        }
        while(true){
            pos=funs.tellg();
            getline(funs,line);
            posTmp=posStruct(line,initFuns.chrIndic,initFuns.zeroBased);
            if(posTmp.chr==posRef.chr){
                if(posTmp.posDNA==posRef.posDNA){
                    if(posTmp.alt.find(posRef.alt)!=std::string::npos){
                        file.at(i)+='\t'+posTmp.score;
                        load(funs,pos);
                        if(posTmp.alt.size()>1){load(funs,pos);}
                        break;
                    }else{
                        if(posTmp.alt.size()==1 && posTmp.alt<posRef.alt){
                            continue;
                        }
                    }
                }else{
                    if(posTmp.posDNA<posRef.posDNA){
                        continue;
                    }
                }
            }
            load(funs,pos);
            file.at(i)+="\tNA";
            break;
        }
    }
    */
