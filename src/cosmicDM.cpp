/*
 * cosmic_DM.cpp
 *
 *  Created on: 3 déc. 2015
 *      Author: drubay
 */


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>	//for sort function
#include <vector>




using namespace std;


int main(int argc, char* argv[]){
	/* import */
	ifstream data(argv[1], std::ios::in);
	ofstream expor(argv[2], ios::out | ios::trunc);
	/*checks*/
	if (!data){cerr << "Input does not exist." << std::endl;return 1;}
	if(!expor){cerr << "Output cannot be created." << endl;return 1;}

	/*Data management*/

        string ldata,elts,chr,ref,alt,posi,sample,conf;
        string saveChr,saveRef,saveAlt,saveSample,savePosi;
        int count=0,sep1,sep2,less,testRec,found;
        int saveTotalRec=1,saveConfRec,saveUnknRec;



	getline(data,ldata,'\n');
	istringstream line(ldata);
	count=0;
	while(getline(line,elts,'\t')){
            switch (count) {
            case 0:
                sep1=elts.find(':');
                sep2=elts.find('-');
                saveChr=elts.substr(0,sep1);
                savePosi=elts.substr(sep1+1,sep2-sep1-1); //as 1-based, same pos at begin-end for SNPs
                break;
            case 1:
                saveRef=elts;
                if(elts.length()!=1){less++;}
                break;
            case 2:
                saveAlt=elts;
                if(elts.length()!=1){less++;}
                break;
            case 3:
                found=elts.find("TCGA");
                if(found>=0){saveSample="TCGA";}else{saveSample="other";} //if missing, -1
                break;
            case 4:
                if(elts=="Confirmed somatic variant"){
                    conf="conf";
                    saveConfRec=1;
                    saveUnknRec=0;
                }else{
                    conf="unkn";
                    saveConfRec=0;
                    saveUnknRec=1;
                }
                break;
            }
            count++;
	}



	while(getline(data,ldata,'\n')){ 
		if(ldata.length()!=0){
			istringstream line(ldata);
			count=0;
			less=0;
                        testRec=0;
			while(getline(line,elts,'\t') && less==0){//elements of chain
                            switch (count) {
                            case 0:
                                sep1=elts.find(':');
                                sep2=elts.find('-');
                                chr=elts.substr(0,sep1);
                                posi=elts.substr(sep1+1,sep2-sep1-1); //as 1-based, same pos at begin-end for SNPs
                                break;
                            case 1:
                                ref=elts;
                                if(elts.length()!=1){less++;}
                                break;
                            case 2:
                                alt=elts;
                                if(elts.length()!=1){less++;}
                                break;
                            case 3:
                                found=elts.find("TCGA");
                                if(found>=0){sample="TCGA";}else{sample="other";} //if missing, -1
                                break;
                            case 4:
                                if(elts=="Confirmed somatic variant"){
                                    conf="conf";
                                }else{
                                    conf="unkn";
                                }
                                break;
                            }
				count++;
			}
			if(less==0){
                                if(chr==saveChr && posi==savePosi && ref==saveRef && alt==saveAlt){
                                    testRec++;
                                    saveTotalRec++;
                                    if(conf=="conf"){
                                       saveConfRec++;
                                    }else{
                                        saveUnknRec++;
                                    }
                                }
                                if(testRec>=1 && sample=="TCGA" && saveSample!="TCGA"){saveSample="TCGA";}
                                if(testRec==0){ // if actual iteration is not recurrence
                                        expor << saveChr << '\t' << savePosi <<'\t'<< saveRef << '\t' << saveAlt << '\t' << saveSample << '\t' << saveConfRec << '\t' << saveUnknRec << '\t' << saveTotalRec; // << savereseq << '\t'  << '\t' << saveTotalRec+1
					if(!data.eof()){expor << endl;}
					//reinitialize and save info of this line which is a new variant (because testrec==0 -> end of recurrence)
                                        saveChr=chr;
                                        savePosi=posi;
                                        saveRef=ref;
                                        saveAlt=alt;
                                        saveSample=sample;
                                        if(conf=="conf"){
                                           saveConfRec=1;
                                           saveUnknRec=0;
                                        }else{
                                            saveConfRec=0;
                                            saveUnknRec=1;
                                        }
                                        saveTotalRec=1;
				}
			}
		}
	}
        if(testRec==0){
            expor << chr << '\t' << posi <<'\t'<< ref << '\t' << alt << '\t' << sample << '\t' << saveConfRec << '\t' << saveUnknRec << '\t' << saveTotalRec << endl;
        }else{
            expor << saveChr << '\t' << savePosi <<'\t'<< saveRef << '\t' << saveAlt << '\t' << saveSample << '\t' << saveConfRec << '\t' << saveUnknRec << '\t' << saveTotalRec << endl;
        }
	return 0;
}


