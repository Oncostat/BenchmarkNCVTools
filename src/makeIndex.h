#ifndef MAKEINDEX_H
#define MAKEINDEX_H
#include "vsfunctions.h"
#include <vector>


vector<streampos> indexChr(fstream &file,struct fileInfos tmpInit,streampos eof){
    vector<streampos> out(24,eof);
    streampos boundChr=1000000000; //small enough to avoid to miss a chr (maybe long...)
    string tmpChr;
    streampos diffPos,nextPos,pos;
    streampos begin=tmpInit.begin,end=tmpInit.end;
    int chr;
    out.at(0)=tmpInit.begin;
    file.clear();
    file.seekg(tmpInit.begin);
    vector<int> vecChr(24);
    if(tmpInit.charSorted){
        //vecChr={0,11,15,16,17,18,19,20,21,1,2,3,4,5,6,7,8,9,10,12,13,14,22,23};
        vecChr={0,9,10,11,12,13,14,15,16,17,18,1,19,20,21,2,3,4,5,6,7,8,22,23};
    }else{
        vecChr={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
    }
    for(int i=1;i<24;++i){
        while(true){
            //if(testEof){begin=tmpInit.begin;testEof=false;}
            nextPos=begin+boundChr;
            if((nextPos+(streampos)15)>eof){
//                if(loopEof){
                    out.at(vecChr[i])=eof;
//                    loopEof=false;
                    break;
//                }else{
//                    testEof=true;
//                    loopEof=true;
//                    continue;
//                }
            }
            file.clear();
            file.seekg(nextPos);
            file.ignore(numeric_limits<streamsize>::max(),'\n');
            pos=file.tellg();
            file>>tmpChr;
            if(tmpInit.chrIndic){tmpChr=tmpChr.substr(3,tmpChr.size());}
            chr=convChr(tmpChr);
            if(chr==vecChr[i]){ //found //begin stay begin
                end=nextPos;
                while(true){ //dichot
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
                out.at(vecChr[i])=pos; //////////////////////////if chr ordered?
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
//    cout << curFilePos << '\t' << tmpChr << '\t' << elt << endl;
    if(curGenPos<stoull(elt)){
        file.clear();
        file.seekg(nextPos);
        return curFilePos;
    }
    while(true){
//        cout << boundTmp<<endl;
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
                if(boundTmp<15){ //if not found because NA -> nextpos is the next, after NAs
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





/*vector<vecPos> indexing(fstream &cadd,fstream &dann,fstream &fath,fstream &funs,struct fileInfos tmpInit,vector<streampos> chrLims,long long int genPos,streampos eof){
    //make vector of genPos for index
    //make an empty struct of vectors of genPos/filePos
    //dichot to find each genPos, if find, report in the empty struct
    vector<long long int> chrSize{249250622,243199374,198022431,191154277,180915261,171115068,159138664,146364023,141213432,135534748,135006517,133851896,115169879,107349541,102531393,90354754,81195211,78077249,59128984,63025521,48129896,51304567,155270561,59373567};
    streampos boundTmp,filePosBound=10000000,
    long long int genPosBound=1000000;
    vector<vector<streamPos>> filePos(4);
    vector<long long int> genPos,curGenPos;
    struct vecPos tmpPos;
    string tmpChr,tmpElt,elt;
    int chr;
    for(int i=0;i<24;++i){
        curGenPos=0;
        while(curGenPos<chrSize.at(i)){
            curGenPos+=genPosBound;
            genPos.push_back(curGenPos);
            filePos.at(0).push_back(getFilePos(cadd))

        }
        tmpPos.filePos.push_back(filePos);
        tmpPos.genPos.push_back(genPos);
    }
}*/






/*
vector<vecPos> indexing(fstream &file,struct fileInfos tmpInit){
    streampos boundChr=1000000000; //small enough to avoid to miss a chr (maybe long...)
    streampos boundPos=10000000;
    streampos bound;
    file.seekg(0, ios::end); //seeking to the end of the file
    streampos eof=file.tellg(); //Getting the length of file or positon of enf of file
    long long int stepGenPos=1000000;
    vector<long long int> chrSize{249250622,243199374,198022431,191154277,180915261,171115068,159138664,146364023,141213432,135534748,135006517,133851896,115169879,107349541,102531393,90354754,81195211,78077249,59128984,63025521,48129896,51304567,155270561,59373567};
    struct vecPos vecTmp;
    vector<vecPos> out;
    string tmpElt,elt,lfile,tmpChr;
    streampos diffPos,nextPos,pos;
    streampos begin=tmpInit.begin,end=tmpInit.end;
    int chr;
    long long int start,nextGenPos;
    bool testNotEof=true;
    for(int i=0;i<24;++i){
        cout << "chr" <<i+1<<'\t';
        begin=tmpInit.begin;
        if(i==0){
            vecTmp.filePos.push_back(tmpInit.begin);
            file.clear();
            file.seekg(tmpInit.begin);
            if(tmpInit.zeroBased){file>>tmpChr>>tmpElt>>elt;}else{file>>tmpChr>>elt;}
            vecTmp.genPos.push_back(stoull(elt));
        }else{
            while(true){
                nextPos=begin+boundChr;
                if((nextPos+(streampos)15)>eof){testNotEof=false;break;}
                file.clear();
                file.seekg(nextPos);
                getline(file,lfile);
                file>>tmpChr;
                if(tmpInit.chrIndic){tmpChr=tmpChr.substr(3,tmpChr.size());}
                chr=convChr(tmpChr);
                getline(file,lfile);
                if(chr==i){ //found //begin stay begin
                    end=nextPos;
                    while(true){ //dichot
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
                        getline(file,lfile);
                        file>>tmpChr;
                        if(tmpInit.chrIndic){tmpChr=tmpChr.substr(3,tmpChr.size());}
                        chr=convChr(tmpChr);
                        if(chr==i){
                            end=nextPos;
                        }else{
                            begin=nextPos;
                        }
                    }
                    vecTmp.filePos.push_back(pos);
                    file.clear();
                    file.seekg(pos);
                    getline(file,lfile);
                    if(tmpInit.zeroBased){file>>tmpChr>>tmpElt>>elt;}else{file>>tmpChr>>elt;}
                    vecTmp.genPos.push_back(stoull(elt));
                    break;
                }else{
                    begin=nextPos;
                }
            }
        }
        //chr found, next pos until the end of this chr
        if(testNotEof){
            begin=vecTmp.filePos[0];
            start=vecTmp.genPos[0];
            nextGenPos=start+stepGenPos;
            while(nextGenPos<chrSize.at(i)){
                bound=boundPos;
                while(true){
                    nextPos=begin+bound;
                    if(nextPos>eof){
                        if(bound>15){
                            bound=(streampos)round(bound/2);
                            continue;
                        }else{
                            vecTmp.filePos.push_back(eof);
                            vecTmp.filePos.push_back(eof);
                            vecTmp.genPos.push_back((long long int)1);
                            vecTmp.genPos.push_back((long long int)chrSize.at(i));
                            break;
                        }
                    }
                    file.clear();
                    file.seekg(nextPos);
                    getline(file,lfile);
                    pos=file.tellg();
                    if(tmpInit.zeroBased){file>>tmpChr>>tmpElt>>elt;}else{file>>tmpChr>>elt;}
                    if(tmpInit.chrIndic){tmpChr=tmpChr.substr(3,tmpChr.size());}
                    chr=convChr(tmpChr);
                    if(chr!=i){bound=(streampos)round(bound/2);if(bound<15){break;}continue;}
                    if(stoull(elt)<nextGenPos){
                        begin=nextPos;
                    }else{
                        if(stoull(elt)>nextGenPos){
                            //if(bound==305){cout << lfile << endl;}
                            //cout << bound << endl;
                            bound=(streampos)round(bound/2);
                            if(bound<15){ //if not found because NA -> nextpos is the next, after NAs
//                                file.clear();
//                                file.seekg(nextPos);
//                                getline(file,lfile);
//                                pos=file.tellg();
                                break;
//                                nextGenPos=stoull(elt);
//                                begin-=10000; //go back to redo loop to find this new position
//                                bound=10000;
                            }
                        }else{
                            end=nextPos;
                            while(true){
                                diffPos=(streampos)(round((end-begin)/2));
                                if(diffPos<15){
                                    file.clear();
                                    file.seekg(begin);
                                    file.ignore(numeric_limits<streamsize>::max(),'\n');
                                    file.ignore(numeric_limits<streamsize>::max(),'\n');
                                    pos=file.tellg();
                                    break;
                                }
                                nextPos=begin+diffPos;
                                file.clear();
                                file.seekg(nextPos);
                                getline(file,lfile);
                                if(tmpInit.zeroBased){file>>tmpChr>>tmpElt>>elt;}else{file>>tmpChr>>elt;}
                                if(stoull(elt)==nextGenPos){
                                    end=nextPos;
                                }else{
                                    begin=nextPos;
                                }
                            }
                            break;
                        }
                    }
                }
                vecTmp.filePos.push_back(pos);
                vecTmp.genPos.push_back(nextGenPos);
                nextGenPos+=stepGenPos;
            }
        }else{
            vecTmp.filePos.push_back(eof);
            vecTmp.filePos.push_back(eof);
            vecTmp.genPos.push_back((long long int)1);
            vecTmp.genPos.push_back((long long int)chrSize.at(i));
        }
        out.push_back(vecTmp);
        vecTmp.filePos.clear();
        vecTmp.genPos.clear();
    }
    //for the end of file, create a 25th chr
    vecTmp.filePos.push_back(eof);
    vecTmp.genPos.push_back((long long int)1);
    out.push_back(vecTmp);
    vecTmp.filePos.clear();
    vecTmp.genPos.clear();
    cout << endl;
    return out;
}
*/


struct posLim posIndex(fstream &file,int refChr){
    struct posLim tmpPos;
    string tmp,line;
    int chr;
    file >> tmp;
    chr=stoi(tmp);
    if(chr==refChr){
        tmpPos.nextPos=chr;
        file >> tmp;
//        cout << tmp << '\t';
        tmpPos.start=stoull(tmp);
        getline(file,line);
        file >> tmp >> line;
//        cout << tmp << endl;
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



