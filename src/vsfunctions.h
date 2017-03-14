#ifndef VSFUNCTIONS_H
#define VSFUNCTIONS_H

using namespace std;


struct position{
    int chr,posDNA;
    string alt,score;
};

struct dataFeatures{
    streampos eofPos;
    bool charSorted;
};


struct fileInfos{
    streampos begin,end;
    bool charSorted,chrIndic,zeroBased;
};

struct init{
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
    getline(elts,elt,'\t'); 
    getline(elts,elt,'\t');
    tmpPos.alt=elt;
    getline(elts,elt,'\t');
    elt.erase(remove(elt.begin(), elt.end(), '\r'), elt.end()); //remove the special dos end of line character ^M (for DANN)
    tmpPos.score=elt;
    return tmpPos;
}




struct dataFeatures searchEoF(fstream &file,streampos begin,streampos stepSize,bool chrIndic){
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
        getline(file,lfile); 
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
    tmpFeat=searchEoF(data,0,1000000000,tmpInit.chrIndic);
    tmpInit.end=tmpFeat.eofPos;
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


vector<vecPos> loadIndex(fstream &file){ 
    struct vecPos vecTmp;
    vector<vecPos> out;
    vector<streampos> filePos(4);
    string tmpChr,genPos,pos;
    int chr,saveChr=0;
    while(!file.eof()){
        file >> tmpChr >> genPos;
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
            vecTmp.genPos.push_back(vecTmp.genPos.at(vecTmp.genPos.size()-1)*2); 
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
    while(true){
        pos=file.tellg();
        getline(file,line);
        posTmp=posStruct(line,chrIndic,zeroBased);
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


#endif // VSFUNCTIONS_H
