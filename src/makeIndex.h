#ifndef MAKEINDEX_H
#define MAKEINDEX_H
#include "vsfunctions.h"
#include <vector>


vector<streampos> indexChr(fstream &file,struct fileInfos tmpInit,streampos eof){
    vector<streampos> out(24,eof);
    streampos boundChr=1000000000;
    string tmpChr;
    streampos diffPos,nextPos,pos;
    streampos begin=tmpInit.begin,end=tmpInit.end;
    int chr;
    out.at(0)=tmpInit.begin;
    file.clear();
    file.seekg(tmpInit.begin);
    vector<int> vecChr(24);
    if(tmpInit.charSorted){
        vecChr={0,9,10,11,12,13,14,15,16,17,18,1,19,20,21,2,3,4,5,6,7,8,22,23};
    }else{
        vecChr={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
    }
    for(int i=1;i<24;++i){
        while(true){
            nextPos=begin+boundChr;
            if((nextPos+(streampos)15)>eof){
                    out.at(vecChr[i])=eof;
                    break;
            }
            file.clear();
            file.seekg(nextPos);
            file.ignore(numeric_limits<streamsize>::max(),'\n');
            pos=file.tellg();
            file>>tmpChr;
            if(tmpInit.chrIndic){tmpChr=tmpChr.substr(3,tmpChr.size());}
            chr=convChr(tmpChr);
            if(chr==vecChr[i]){
                end=nextPos;
                while(true){
                    diffPos=(streampos)(round((end-begin)/2));
                    nextPos=begin+diffPos;
                    if(diffPos<15){
                        file.clear();
                        file.seekg(begin);
                        file.ignore(numeric_limits<streamsize>::max(),'\n');
                        file.ignore(numeric_limits<streamsize>::max(),'\n');
                        pos=file.tellg();
                        break;
                    }
                    file.clear();
                    file.seekg(nextPos);
                    file.ignore(numeric_limits<streamsize>::max(),'\n');
                    file>>tmpChr;
                    if(tmpInit.chrIndic){tmpChr=tmpChr.substr(3,tmpChr.size());}
                    chr=convChr(tmpChr);
                    if(chr==vecChr[i]){
                        end=nextPos;
                    }else{
                        begin=nextPos;
                    }
                }
                out.at(vecChr[i])=pos;
                break;
            }else{
                begin=nextPos;
            }
        }
    }
    out.push_back(eof);
    return out;
}




streampos getFilePos(long long int curGenPos,fstream &file,bool zeroBased,streampos endChr,streampos eof,streampos curFilePos){
    streampos boundTmp=10000000;
    streampos nextPos,diffPos;
    string tmpChr,tmpElt,elt;
    if(curFilePos==eof){
        return eof;
    }
    file.clear();
    file.seekg(curFilePos);
    if(zeroBased){file>>tmpChr>>tmpElt>>elt;}else{file>>tmpChr>>elt;}
    if(curGenPos<stoull(elt)){
        file.clear();
        file.seekg(nextPos);
        return curFilePos;
    }
    while(true){
        nextPos=curFilePos+boundTmp;
        if(nextPos>eof || nextPos>endChr){
            if(boundTmp>15){
                boundTmp=(streampos)round(boundTmp/2);
                continue;
            }else{
                if(nextPos>eof){
                    return eof;
                }else{
                    return endChr;
                }
            }
        }
        file.clear();
        file.seekg(nextPos);
        file.ignore(numeric_limits<streamsize>::max(),'\n');
        if(zeroBased){file>>tmpChr>>tmpElt>>elt;}else{file>>tmpChr>>elt;}
        if(stoull(elt)<curGenPos){
            curFilePos=nextPos;
        }else{
            if(stoull(elt)>curGenPos){
                boundTmp=(streampos)round(boundTmp/2);
                if(boundTmp<15){
                    file.clear();
                    file.seekg(curFilePos);
                    file.ignore(numeric_limits<streamsize>::max(),'\n');
                    file.ignore(numeric_limits<streamsize>::max(),'\n');
                    return file.tellg();
                }
            }else{
                endChr=nextPos;
                while(true){
                    diffPos=(streampos)(round((endChr-curFilePos)/2));
                    if(diffPos<15){
                        file.clear();
                        file.seekg(curFilePos);
                        file.ignore(numeric_limits<streamsize>::max(),'\n');
                        file.ignore(numeric_limits<streamsize>::max(),'\n');
                        return file.tellg();
                    }
                    nextPos=curFilePos+diffPos;
                    file.clear();
                    file.seekg(nextPos);
                    file.ignore(numeric_limits<streamsize>::max(),'\n');
                    if(zeroBased){file>>tmpChr>>tmpElt>>elt;}else{file>>tmpChr>>elt;}
                    if(stoull(elt)==curGenPos){
                        endChr=nextPos;
                    }else{
                        curFilePos=nextPos;
                    }
                }
            }
        }
    }
}





struct posLim posIndex(fstream &file,int refChr){
    struct posLim tmpPos;
    string tmp,line;
    int chr;
    file >> tmp;
    chr=stoi(tmp);
    if(chr==refChr){
        tmpPos.nextPos=chr;
        file >> tmp;
        tmpPos.start=stoull(tmp);
        getline(file,line);
        file >> tmp >> line;
        tmpPos.stop=stoull(line);
        return tmpPos;
    }else{
        file.clear();
        file.seekg(0);
        while(getline(file,line)){
            file >> tmp;
            chr=stoi(tmp);
            if(chr==refChr){
                tmpPos.nextPos=chr;
                file >> tmp;
                tmpPos.start=stoull(tmp);
                getline(file,line);
                file >> tmp >> line;
                tmpPos.stop=stoull(line);
                return tmpPos;
            }
        }
    }
    cout << "This chromosome does not exist in this database: " << refChr <<endl;
    tmpPos.nextPos=-2;
    return tmpPos;
}



#endif // MAKEINDEX_H



