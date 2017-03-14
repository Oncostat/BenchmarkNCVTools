#ifndef VSFUNCTIONS_H
#define VSFUNCTIONS_H
//#include <iostream>

using namespace std;


struct position{
    int chr,posDNA;
    string alt,score;
//    char search;
/*    position(){
        chr=0;posDNA=0;alt="NA";score="NA";
    }*/
//    position(int a,string b, string c):posDNA(a),alt(b),line(c){}
};

struct dataFeatures{
    streampos eofPos;
    bool charSorted;
//    dataFeatures():filePos(24){}
//    vector<streampos> filePos;
};


struct fileInfos{ //detect the number of columns until the first letter? -> if 2 columns of coordinates (ref is at the 4th column), it is zerobased -> add a bool zeroBased
    streampos begin,end;
    bool charSorted,chrIndic,zeroBased;
};

struct init{ //detect the number of columns until the first letter? -> if 2 columns of coordinates (ref is at the 4th column), it is zerobased -> add a bool zeroBased
    streampos begin,end;
    bool chrIndic,zeroBased;
    vector<streampos> chrPos;
    init(){
        for(int i=0;i<24;++i){
            chrPos.push_back(-1);
        }
    }
};

struct posLim{
    streampos start,stop;
    long long int nextPos;
};



inline int convChr(string chr){
    if(chr=="X"){
        return 22;
    }else{
        if(chr=="Y"){
            return 23;
        }else{
            int tmp;
            tmp=stoi(chr)-1;
            if(tmp<24){
                return tmp;
            }else{
                return 9999;
            }
        }
    }
}



inline struct position posStruct(string line,bool chrIndic,bool zeroBased){
    struct position tmpPos;
    string elt;
    istringstream elts(line);
    getline(elts,elt,'\t');
    if(chrIndic){elt=elt.substr(3,elt.size());}
    tmpPos.chr=convChr(elt);
    getline(elts,elt,'\t');
    if(zeroBased){
        getline(elts,elt,'\t');
        tmpPos.posDNA=stoull(elt);
    }else{
        tmpPos.posDNA=stoull(elt);
    }
    getline(elts,elt,'\t'); //skip ref
    getline(elts,elt,'\t');
    tmpPos.alt=elt;
    getline(elts,elt,'\t');
    elt.erase(remove(elt.begin(), elt.end(), '\r'), elt.end()); //remove the special dos end of line character ^M (for DANN)
    tmpPos.score=elt;
    return tmpPos;
}




struct dataFeatures searchEoF(fstream &file,streampos begin,streampos stepSize,bool chrIndic){ ///modify to search also strange chr which will be the eof too
    struct dataFeatures tmpFeat;
    vector<string> vecChr{"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","X","Y"};
    vector<int> testSort(vecChr.size(),0);
    int rank=1;
    int tmpPos=0;
    vector <string>::iterator found;
    string elt,tmpElt;
    string lfile;
    streampos nextPos;
    bool eofFound=true;
    file.seekg(0, ios::end);
    tmpFeat.eofPos=file.tellg();
    file.seekg(0);
    while(eofFound){
        nextPos=begin+stepSize;
        if(nextPos>tmpFeat.eofPos){break;}
        file.clear();
        file.seekg(nextPos);
        getline(file,lfile); // in the case of pointer in the middle of a line
        getline(file,lfile);
        istringstream elts(lfile);
        getline(elts,elt,'\t');
        if(!file.eof()){
            if(chrIndic){
                tmpElt=elt.substr(3,elt.size());
                found=find(vecChr.begin(),vecChr.end(),tmpElt);
            }else{
                found=find(vecChr.begin(),vecChr.end(),elt);
            }
            if(found!=vecChr.end()){
                tmpPos=distance(vecChr.begin(),found);
                if(testSort[tmpPos]==0){
                    testSort[tmpPos]=rank;
                    ++rank;
                }
            }
        }

        if(file.eof() || found==vecChr.end()){
            if(stepSize<15){
                getline(file,lfile);
                tmpFeat.eofPos=file.tellg();
                eofFound=false;
            }else{
                stepSize=(streampos)round(stepSize/2);
            }
        }else{
            begin=nextPos;
        }
    }
    tmpFeat.charSorted=false;
    rank=0;
    for(int i=0;i<testSort.size();++i){
        if(testSort[i]!=0){
            ++rank;
            if(testSort[i]!=rank){
                tmpFeat.charSorted=true;
                break;
            }
        }
    }
    return tmpFeat;
}




struct fileInfos initializeData(fstream &data,string nameData){
    struct fileInfos tmpInit;
    struct dataFeatures tmpFeat;
    string ldata,elt;
    streampos savePos=data.tellg();
    vector<string>::iterator nuclFound;
    vector<string> nucl{"A","C","G","T"};
    while(getline(data,ldata)){
        if(ldata.find("#")==string::npos){
            istringstream line(ldata);
            getline(line,elt,'\t');
            int chrFound=elt.find("chr");
            getline(line,elt,'\t');
            getline(line,elt,'\t');
            nuclFound=find(nucl.begin(),nucl.end(),elt);
            if(nuclFound!=nucl.end()){
                tmpInit.zeroBased=false;
            }else{
                getline(line,elt,'\t');
                nuclFound=find(nucl.begin(),nucl.end(),elt);
                if(nuclFound!=nucl.end()){
                    tmpInit.zeroBased=true;
                }else{
                    cerr << "The format of " << nameData << " cannot be determined. Please check if you have the format: chr\tposition\treference\talternative\t...\nThe position should be 1 based. For zero based data, the format chr\tposition0\tposition1\treference\talternative\t... is supported" << endl;
                }
            }
            tmpInit.begin=savePos;
            if(chrFound==string::npos){
                tmpInit.chrIndic=false;
            }else{
                tmpInit.chrIndic=true;
            }
            break;
        }
        savePos=data.tellg();
    }
/*    data.seekg(0, ios::end);
    cout << data.tellg()<<endl;
    data.seekg(0);
    cout << data.tellg()<<endl;
    cout << "GO SEARCH\n";*/
    tmpFeat=searchEoF(data,0,1000000000,tmpInit.chrIndic);
    /*cout << "END SEARCH\n";
    data.seekg(0, ios::end);
    cout << data.tellg()<<endl;
    data.seekg(0);
    cout << data.tellg()<<endl;*/
    tmpInit.end=tmpFeat.eofPos; //find the position of the last line
    tmpInit.charSorted=tmpFeat.charSorted;
    return tmpInit;
}


inline void load(fstream &file,streampos pos){
    file.clear();
    file.seekg(pos);
}

struct fileInfos fromInit(fstream &file){
    struct fileInfos tmpInit;
    string fileName,begin,end,chrIndic,zeroBased;
    file >> fileName >> begin >> end >> chrIndic >> zeroBased;
    tmpInit.begin=stoull(begin);
    tmpInit.end=stoull(end);
    tmpInit.chrIndic=(chrIndic=="1");
    tmpInit.zeroBased=(zeroBased=="1");
    return tmpInit;
}


struct vecPos{
    vector<vector<streampos>> filePos;
    vector<long long int> genPos;
};


vector<vecPos> loadIndex(fstream &file){ ///here adapt to new filePos
    //if all posfile equal to previous, don't push back
    struct vecPos vecTmp;
    vector<vecPos> out;
    vector<streampos> filePos(4);
    string tmpChr,genPos,pos;
    int chr,saveChr=0;
    while(!file.eof()){
        file >> tmpChr >> genPos;
//        cout << "test\t"<<tmpChr << '\t' << genPos << '\n';
//        cout.flush();
        for(int i=0;i<4;++i){
            file >> pos;
            filePos.at(i)=stoull(pos);
        }
        chr=stoi(tmpChr);
        if(chr==saveChr){
            vecTmp.filePos.push_back(filePos);
            vecTmp.genPos.push_back(stoull(genPos));
        }else{
            vecTmp.filePos.push_back(filePos);
            vecTmp.genPos.push_back(vecTmp.genPos.at(vecTmp.genPos.size()-1)*2); //add last line of chr which is not limited by posDNA, but by filePos in main.cpp
            out.push_back(vecTmp);
            vecTmp.filePos.clear();
            vecTmp.genPos.clear();
            vecTmp.filePos.push_back(filePos);
            vecTmp.genPos.push_back(stoull(genPos));
            saveChr=chr;
        }
    }
    return out;
}

inline void scoring(struct position * posRef, fstream &file,bool chrIndic,bool zeroBased,string* toFill){
    struct position posTmp;
    string line;
    streampos pos=file.tellg();
//    if(posRef->chr==0 && posRef->posDNA==120302485){cout << posRef->chr << ":" << posRef->posDNA << "-" << posRef->alt<<endl;cout << pos << endl;}
    while(true){
        pos=file.tellg();
        getline(file,line);
//      if(posRef->chr==0 && posRef->posDNA==120302485){cout <<line <<endl;}
        posTmp=posStruct(line,chrIndic,zeroBased);
//        if(posRef->chr==0 && posRef->posDNA==120302485){ cout << posTmp.chr << ":" << posTmp.posDNA << "-" << posTmp.alt<<endl;}//ref=A alt=T
        if(posTmp.chr==posRef->chr){
            if(posTmp.posDNA==posRef->posDNA){
                if(posTmp.alt.find(posRef->alt)!=std::string::npos){
                    *toFill+='\t'+posTmp.score;
                    load(file,pos);
                    return;
                }else{
                    if(posTmp.alt.size()==1 && posTmp.alt<posRef->alt){
                        continue;
                    }
                }
            }else{
                if(posTmp.posDNA<posRef->posDNA){
                    continue;
                }
            }
        }
        load(file,pos);
        *toFill+="\tNA";
        return;
    }
}

/*vector<vecPos> loadIndex(fstream &file){ ///here adapt to new filePos
    //if all posfile equal to previous, don't push back
    struct vecPos vecTmp;
    vector<vecPos> out;
    dataInIndex toFill;
    string tmpChr,genPos,filePos;
    int chr,saveChr=0;
    while(!file.eof()){
        file >> tmpChr >> genPos >> filePos;
        chr=stoi(tmpChr);
        if(chr==saveChr){
            vecTmp.filePos.push_back(stoull(filePos));
            vecTmp.genPos.push_back(stoull(genPos));
            vecTmp.classified.push_back(toFill);
        }else{
            vecTmp.filePos.push_back(stoull(filePos));
            vecTmp.genPos.push_back(vecTmp.genPos.at(vecTmp.genPos.size()-1)*2); //add last line of chr which is not limited by posDNA, but by filePos in main.cpp
            vecTmp.classified.push_back(toFill);
            out.push_back(vecTmp);
            vecTmp.filePos.clear();
            vecTmp.genPos.clear();
            vecTmp.filePos.push_back(stoull(filePos));
            vecTmp.genPos.push_back(stoull(genPos));
            vecTmp.classified.push_back(toFill);
            saveChr=chr;
        }
    }
    return out;
}*/


/*
void updateThres(struct posLim * lim,int fileNumber,vector<vecPos> * pos,int posVecChr,int posVecDNA){
    lim->start=pos->at(posVecChr).filePos[posVecDNA].at(fileNumber);
    lim->stop=pos->at(posVecChr).filePos[posVecDNA+1].at(fileNumber);
    lim->nextPos=pos->at(posVecChr).genPos[posVecDNA+1].at(fileNumber);
}*/



/*
char tryPos(struct position *ref,string line,bool chrIndic,bool zeroBased){
    if(tmpPos->chr==9999){
        return 'n';//no output
    }
    struct position tmpPos=posStruct(line,chrIndic,zeroBased);
    if(tmpPos->chr==ref->chr){
        if(tmpPos->posDNA==ref->posDNA){
            if(tmpPos->alt.find(ref->alt)!=std::string::npos){
                return 'h';
            }else{
                return 'c';
            }
        }else{
            return 'n';
        }
    }else{
        return 'n';
    }
}
*/




/*struct posLim posFile(fstream &file,long long int refPos,struct posLim posThres,string nameFile){
    struct posLim tmpPos;
    string tmp,line;
    long long int position;
    streampos pos=posThres.start;
    if(file.tellg()<pos && file.tellg()>=posThres.stop){
        file.clear();
        file.seekg(pos);
    }else{
        pos=file.tellg();
    }
    file >> tmp;
    position=stoull(tmp);
//    cout << refPos << '\t' << position << '\t' << pos << '\t'<< posThres.stop <<endl;
    if(position==refPos){
        tmpPos.pos=position;
        file >> tmp;
        tmpPos.start=stoull(tmp);
        getline(file,line);
        file >> tmp >> line;
        tmpPos.stop=stoull(line);
        return tmpPos;
    }else{
        while(getline(file,line) && pos<posThres.stop){
            pos=file.tellg();
            file >> tmp;
            position=stoull(tmp);
            if(position>refPos){
                file.clear();
                file.seekg(pos);
                file >> tmp;
                position=stoull(tmp);
                tmpPos.pos=position;
                file >> tmp;
                tmpPos.start=stoull(tmp);
                getline(file,line);
                file >> tmp >> line;
                tmpPos.stop=stoull(line);
                return tmpPos;
            }
        }
    }
//    cout << posThres.start << '\t' << file.tellg() << '\t' << posThres.stop << '\t' << refPos << '\t' << position << endl;
    cout << "This position does not exist in " << nameFile << " database: chromosome " << posThres.pos+1 << " position "<< refPos <<endl;
    tmpPos.pos=-2;
    return tmpPos;
}*/


#endif // VSFUNCTIONS_H
