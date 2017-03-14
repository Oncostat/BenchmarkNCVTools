/*
 * mergesnpsom.cpp
 *
 *  Created on: 5 janv. 2016
 *      Author: drubay
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <map>



using namespace std;

int main(int argc, char* argv[])
{
	ifstream snp(argv[1], std::ios::in);
	ifstream som(argv[2], std::ios::in);
	ofstream output(argv[3], ios::out | ios::trunc);

	if(!snp){cerr << "The SNP file does not exist." << endl;return 1;}
	if(!som){cerr << "The SOM file does not exist." << endl;return 1;}
	if(!output){cerr << "The output file cannot be created." << std::endl;return 1;}

	int count,snpchr,snpbegin,snpend,somchr,sombegin,somend,test,begin;
	string cllsc,liversc,lungsc,melsc;
	string lsnp,lsom,elts,snpsc;
	streampos pos=0;

	map< string , int > chrMap;
	map< int , string > mapChr;

	for(int i=0;i<25;i++){
		if(i<23){
			chrMap[to_string(i)]=i;
			mapChr[i]=to_string(i);
		}else{
			if(i==23){
				chrMap["X"]=i;
				mapChr[i]="X";
			}else{
				chrMap["Y"]=i;
				mapChr[i]="Y";
			}
		}
	}


	getline(snp,lsnp,'\n');
	istringstream eltssnp(lsnp);
	count=0;
	while(getline(eltssnp,elts,'\t')){
		if(count==0){snpchr=chrMap[elts];} //istringstream(elts) >> snpchr;
		if(count==1){istringstream(elts) >> snpbegin;}
		if(count==2){istringstream(elts) >> snpend;}
		if(count==3){snpsc=elts;}
		count++;
	}

	getline(som,lsom,'\n');
	istringstream eltssom(lsom);
	count=0;
	while(getline(eltssom,elts,'\t')){
		if(count==0){somchr=chrMap[elts];} //istringstream(elts) >> somchr;}
		if(count==1){istringstream(elts) >> sombegin;}
		if(count==2){istringstream(elts) >> somend;}
		if(count==3){cllsc=elts;}
		if(count==4){liversc=elts;}
		if(count==5){lungsc=elts;}
		if(count==6){melsc=elts;}
		count++;
	}
	som.clear();som.seekg(pos);
        output << lsnp << '\t' << cllsc << '\t' << liversc << '\t' << lungsc << '\t' << melsc << endl; //1st end snp < 1st end som

	while(getline(snp,lsnp,'\n')){
		if(som.eof()){
			output << lsnp << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << endl;
		}else{
			istringstream eltssnp(lsnp);
			count=0;
			while(getline(eltssnp,elts,'\t')){
				if(count==0){snpchr=chrMap[elts];} //istringstream(elts) >> snpchr;}
				if(count==1){istringstream(elts) >> snpbegin;}
				if(count==2){istringstream(elts) >> snpend;}
				if(count==3){snpsc=elts;}
				count++;
			}
			test=0;
			while(test==0){
				pos = som.tellg();
				getline(som,lsom,'\n');
				istringstream eltssom(lsom);
				count=0;
				while(getline(eltssom,elts,'\t')){
					if(count==0){somchr=chrMap[elts];} //istringstream(elts) >> somchr;}
					if(count==1){istringstream(elts) >> sombegin;}
					if(count==2){istringstream(elts) >> somend;}
					if(count==3){cllsc=elts;}
					if(count==4){liversc=elts;}
					if(count==5){lungsc=elts;}
					if(count==6){melsc=elts;}
					count++;
				}
				if(somchr>snpchr){
					output << lsnp << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << endl;
					test++;
					if(!som.eof()){som.clear();som.seekg(pos);}
				}else{
					if(somchr==snpchr){
						if(sombegin>snpend){
							output << lsnp << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << endl;
							test++;
							if(!som.eof()){som.clear();som.seekg(pos);}
						}else{
							if(snpbegin>sombegin){begin=snpbegin;}else{begin=sombegin;}
							if(somend>snpend){
								output << mapChr[snpchr] << '\t' << begin << '\t' << snpend << '\t' << snpsc << '\t' << cllsc << '\t' << liversc << '\t' << lungsc << '\t' << melsc << endl;
								test++;
								som.clear();som.seekg(pos);
							}else{
								if(somend==snpbegin){
									output << mapChr[snpchr] << '\t' << somend << '\t' << somend+1 << '\t' << snpsc << '\t' << cllsc << '\t' << liversc << '\t' << lungsc << '\t' << melsc << endl;
								}else{
									if(somend==snpend){
										output << mapChr[snpchr] << '\t' << begin << '\t' << somend << '\t' << snpsc << '\t' << cllsc << '\t' << liversc << '\t' << lungsc << '\t' << melsc << endl; // << setprecision(20)
										test++;
									}else{
										if(somend>snpbegin && somend<snpend){
											output << mapChr[snpchr] << '\t' << begin << '\t' << somend << '\t' << snpsc << '\t' << cllsc << '\t' << liversc << '\t' << lungsc << '\t' << melsc << endl; // << setprecision(20)
										}
									}
								}
							}
						}
					}
				}
				if(som.eof()){test++;}
			}
		}
	}
	snp.close();
	som.close();
	return 0;
}


